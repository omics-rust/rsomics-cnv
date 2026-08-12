use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::signals::{RequiredSignals, SampleSelection, SignalReader};

pub const BINS: usize = 150;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolysomyOptions {
    pub fit_threshold: f64,
    pub copy_number_penalty: f64,
    pub peak_symmetry: f64,
    pub minimum_peak_size: f64,
    pub minimum_fraction: f64,
    pub include_aa: bool,
    pub smoothing: i32,
}

impl Default for PolysomyOptions {
    fn default() -> Self {
        Self {
            fit_threshold: 3.3,
            copy_number_penalty: 0.7,
            peak_symmetry: 0.5,
            minimum_peak_size: 0.1,
            minimum_fraction: 0.1,
            include_aa: false,
            smoothing: -3,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct DistributionResult {
    pub sample: String,
    pub chromosomes: Vec<ChromosomeDistribution>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ChromosomeDistribution {
    pub reference_name: String,
    pub observations: u64,
    pub bins: Vec<DistributionBin>,
    pub rr_boundary: usize,
    pub heterozygous_center: usize,
    pub aa_boundary: usize,
    pub fitted_end: usize,
    pub preliminary_copy_number: Option<i8>,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct DistributionBin {
    pub baf: f64,
    pub normalized_count: f64,
}

struct Histogram {
    reference_name: String,
    observations: u64,
    counts: [f64; BINS],
}

impl Histogram {
    fn new(reference_name: String) -> Self {
        Self {
            reference_name,
            observations: 0,
            counts: [0.0; BINS],
        }
    }

    fn push(&mut self, baf: f64) {
        let index = (baf * (BINS - 1) as f64) as usize;
        self.counts[index.min(BINS - 1)] += 1.0;
        self.observations += 1;
    }
}

pub fn analyze_distributions(
    input: &Path,
    sample: Option<String>,
    options: PolysomyOptions,
) -> Result<DistributionResult> {
    validate_options(options)?;
    let mut reader = SignalReader::open(
        input,
        SampleSelection {
            query: sample,
            control: None,
        },
        RequiredSignals::Baf,
    )?;
    let sample = reader.query_sample().to_owned();
    let mut histograms = Vec::new();
    while let Some(site) = reader.next_site()? {
        if histograms
            .last()
            .is_none_or(|histogram: &Histogram| histogram.reference_name != site.reference_name)
        {
            histograms.push(Histogram::new(site.reference_name));
        }
        let baf = site
            .query
            .baf
            .ok_or_else(|| invalid("signal reader returned a site without query BAF"))?;
        histograms.last_mut().unwrap().push(baf);
    }
    if histograms.is_empty() {
        return Err(invalid(
            "no informative BAF sites remain after sample selection",
        ));
    }
    let chromosomes = histograms
        .into_iter()
        .map(|histogram| normalize(histogram, options.smoothing))
        .collect::<Result<Vec<_>>>()?;
    Ok(DistributionResult {
        sample,
        chromosomes,
    })
}

fn normalize(histogram: Histogram, smoothing: i32) -> Result<ChromosomeDistribution> {
    let window = smoothing_window(smoothing)?;
    let smoothed = smooth(&histogram.counts, window);
    let mut rr_boundary = 0;
    for index in 0..BINS / 2 {
        if smoothed[index] < smoothed[rr_boundary] {
            rr_boundary = index;
        }
    }
    let mut aa_boundary = BINS - 1;
    for index in (BINS / 2..BINS).rev() {
        if smoothed[index] < smoothed[aa_boundary] {
            aa_boundary = index;
        }
    }
    rr_boundary += window / 2;
    aa_boundary = (aa_boundary + window / 2).min(BINS - 1);
    if rr_boundary >= aa_boundary {
        return Err(invalid(format!(
            "distribution normalization failed for {}: boundaries {rr_boundary} and {aa_boundary}",
            histogram.reference_name
        )));
    }

    let mut counts = if smoothing > 0 {
        smoothed
    } else {
        histogram.counts.to_vec()
    };
    let mut fitted_end = aa_boundary;
    for index in aa_boundary..BINS {
        if counts[fitted_end] < counts[index] {
            fitted_end = index;
        }
    }

    let (mut rr_max, mut ra_max, mut aa_max) = (0.0f64, 0.0f64, 0.0f64);
    let (mut rr_sum, mut ra_sum, mut aa_sum) = (0.0f64, 0.0f64, 0.0f64);
    for value in &counts[..rr_boundary] {
        rr_sum += value;
        rr_max = rr_max.max(*value);
    }
    for value in &counts[rr_boundary..=aa_boundary] {
        ra_sum += value;
        ra_max = ra_max.max(*value);
    }
    for value in &counts[aa_boundary + 1..] {
        aa_sum += value;
        aa_max = aa_max.max(*value);
    }
    let preliminary_copy_number =
        if ra_sum == 0.0 || (ra_sum / rr_sum < 0.1 && aa_sum / ra_sum > 1.0) {
            ra_max = aa_max;
            Some(1)
        } else if ra_sum / rr_sum < 0.1 || aa_sum / ra_sum > 1.0 {
            ra_max = aa_max;
            Some(-1)
        } else {
            None
        };
    normalize_range(&mut counts[..rr_boundary], rr_max);
    normalize_range(&mut counts[rr_boundary..=aa_boundary], ra_max);
    normalize_range(&mut counts[aa_boundary + 1..], aa_max);

    Ok(ChromosomeDistribution {
        reference_name: histogram.reference_name,
        observations: histogram.observations,
        bins: counts
            .into_iter()
            .enumerate()
            .map(|(index, normalized_count)| DistributionBin {
                baf: index as f64 / (BINS - 1) as f64,
                normalized_count,
            })
            .collect(),
        rr_boundary,
        heterozygous_center: BINS / 2,
        aa_boundary,
        fitted_end,
        preliminary_copy_number,
    })
}

fn smooth(counts: &[f64; BINS], window: usize) -> Vec<f64> {
    let mut output = vec![0.0; BINS];
    let mut half = window / 2;
    let mut average = counts[0];
    output[0] = average;
    for index in 1..half {
        average += counts[2 * index - 1];
        output[index] = average / (2 * index + 1) as f64;
    }
    average = 0.0;
    for index in 0..BINS {
        average += counts[index];
        if index >= window - 1 {
            output[index - half] = average / window as f64;
            average -= counts[index + 1 - window];
        }
    }
    for index in BINS - half..BINS {
        average -= counts[index - half];
        half -= 1;
        output[index] = average / (2 * half + 1) as f64;
        average -= counts[index - half];
    }
    output
}

fn normalize_range(values: &mut [f64], maximum: f64) {
    if maximum != 0.0 {
        for value in values {
            *value /= maximum;
        }
    }
}

fn validate_options(options: PolysomyOptions) -> Result<()> {
    for (name, value) in [
        ("copy-number penalty", options.copy_number_penalty),
        ("peak symmetry", options.peak_symmetry),
        ("minimum peak size", options.minimum_peak_size),
        ("minimum fraction", options.minimum_fraction),
    ] {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(invalid(format!(
                "{name} must be finite and between 0 and 1"
            )));
        }
    }
    if !options.fit_threshold.is_finite() || options.fit_threshold <= 0.0 {
        return Err(invalid("fit threshold must be finite and positive"));
    }
    smoothing_window(options.smoothing)?;
    Ok(())
}

fn smoothing_window(smoothing: i32) -> Result<usize> {
    let half = if smoothing == 0 {
        3
    } else {
        usize::try_from(smoothing.unsigned_abs())
            .map_err(|_| invalid("smoothing half-width exceeds usize"))?
    };
    if half >= BINS / 2 {
        return Err(invalid(format!(
            "smoothing half-width must be smaller than {}",
            BINS / 2
        )));
    }
    Ok(half * 2 + 1)
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
