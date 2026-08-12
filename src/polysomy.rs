use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::fitting::{FitOutcome, Peak, fit};
use crate::selection::{CompiledSelection, SiteSelection};
use crate::signals::{RequiredSignals, SampleSelection, SelectedSignalReader};

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

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct PolysomyResult {
    pub sample: String,
    pub chromosomes: Vec<ChromosomePolysomy>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ChromosomePolysomy {
    pub distribution: ChromosomeDistribution,
    pub copy_number: f64,
    pub absolute_deviation: Option<f64>,
    pub candidates: Vec<CandidateFit>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct CandidateFit {
    pub model_copy_number: u8,
    pub estimated_copy_number: f64,
    pub absolute_deviation: Option<f64>,
    pub selected: bool,
    pub rejection: Option<FitRejection>,
    pub curves: Vec<FitCurve>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum FitRejection {
    FitThreshold,
    PeakSymmetry,
    PeakSize,
    PeakPlacement,
    CopyNumberPenalty,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct FitCurve {
    pub start_bin: usize,
    pub end_bin: usize,
    pub absolute_deviation: f64,
    pub function: String,
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
        let index = (baf as f32 * (BINS - 1) as f32) as usize;
        self.counts[index.min(BINS - 1)] += 1.0;
        self.observations += 1;
    }
}

pub fn analyze(
    input: &Path,
    sample: Option<String>,
    options: PolysomyOptions,
) -> Result<PolysomyResult> {
    analyze_selected(input, sample, options, &SiteSelection::default())
}

/// Fits chromosome copy-number models after applying regions or targets.
pub fn analyze_selected(
    input: &Path,
    sample: Option<String>,
    options: PolysomyOptions,
    selection: &SiteSelection,
) -> Result<PolysomyResult> {
    let distributions = analyze_distributions_selected(input, sample, options, selection)?;
    let chromosomes = distributions
        .chromosomes
        .into_iter()
        .map(|distribution| fit_distribution(distribution, options))
        .collect::<Result<Vec<_>>>()?;
    Ok(PolysomyResult {
        sample: distributions.sample,
        chromosomes,
    })
}

pub fn analyze_distributions(
    input: &Path,
    sample: Option<String>,
    options: PolysomyOptions,
) -> Result<DistributionResult> {
    analyze_distributions_selected(input, sample, options, &SiteSelection::default())
}

/// Builds chromosome BAF distributions after applying regions or targets.
pub fn analyze_distributions_selected(
    input: &Path,
    sample: Option<String>,
    options: PolysomyOptions,
    selection: &SiteSelection,
) -> Result<DistributionResult> {
    validate_options(options)?;
    let selection = CompiledSelection::new(selection)?;
    let mut reader = SelectedSignalReader::open(
        input,
        SampleSelection {
            query: sample,
            control: None,
        },
        RequiredSignals::Baf,
        selection,
    )?;
    let sample = reader.query_sample().to_owned();
    let mut histograms = Vec::new();
    reader.visit(false, |site| {
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
        Ok(())
    })?;
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

fn fit_distribution(
    distribution: ChromosomeDistribution,
    options: PolysomyOptions,
) -> Result<ChromosomePolysomy> {
    if let Some(copy_number) = distribution.preliminary_copy_number {
        return Ok(ChromosomePolysomy {
            distribution,
            copy_number: f64::from(copy_number),
            absolute_deviation: None,
            candidates: Vec::new(),
        });
    }
    let cn2 = fit_cn2(&distribution, options)?;
    let cn3 = fit_cn3(&distribution, options, cn2.aa.as_ref())?;
    let cn4 = fit_cn4(&distribution, options)?;
    let mut candidates = vec![cn2.candidate, cn3, cn4];
    let mut selected = None;
    let mut best_deviation = candidates[0].absolute_deviation.unwrap_or(f64::INFINITY);
    if candidates[0].rejection.is_none() {
        selected = Some(0);
    }
    for index in 1..candidates.len() {
        if candidates[index].rejection.is_some() {
            continue;
        }
        let deviation = candidates[index]
            .absolute_deviation
            .unwrap_or(f64::INFINITY);
        if selected.is_none() || deviation < (1.0 - options.copy_number_penalty) * best_deviation {
            if let Some(previous) = selected.replace(index) {
                candidates[previous].rejection = Some(FitRejection::CopyNumberPenalty);
            }
            best_deviation = deviation;
        } else {
            candidates[index].rejection = Some(FitRejection::CopyNumberPenalty);
        }
    }
    let (copy_number, absolute_deviation) = if let Some(index) = selected {
        candidates[index].selected = true;
        (
            candidates[index].estimated_copy_number,
            candidates[index].absolute_deviation,
        )
    } else {
        (-1.0, candidates[0].absolute_deviation)
    };
    Ok(ChromosomePolysomy {
        distribution,
        copy_number,
        absolute_deviation,
        candidates,
    })
}

struct Cn2Fit {
    candidate: CandidateFit,
    aa: Option<FitCurve>,
}

fn fit_cn2(distribution: &ChromosomeDistribution, options: PolysomyOptions) -> Result<Cn2Fit> {
    let x = x_values(distribution);
    let y = y_values(distribution);
    let aa = if options.include_aa {
        let peak = Peak::exponential(1.0, 1.0, 0.2, 5)
            .scan(2, 0.01, 0.3)?
            .scan(0, 0.05, 1.0)?;
        let outcome = fit(
            vec![peak],
            &x[distribution.aa_boundary..=distribution.fitted_end],
            &y[distribution.aa_boundary..=distribution.fitted_end],
        )?;
        Some(curve(
            outcome,
            distribution.aa_boundary,
            distribution.fitted_end,
        ))
    } else {
        None
    };
    let peak = Peak::bounded_gaussian(1.0, 0.5, 0.03, 0.45, 0.55, 7)?
        .scan(2, 0.01, 0.3)?
        .scan(0, 0.05, 1.0)?;
    let outcome = fit(
        vec![peak],
        &x[distribution.rr_boundary..=distribution.aa_boundary],
        &y[distribution.rr_boundary..=distribution.aa_boundary],
    )?;
    let ra = curve(outcome, distribution.rr_boundary, distribution.aa_boundary);
    let deviation =
        ra.absolute_deviation + aa.as_ref().map_or(0.0, |curve| curve.absolute_deviation);
    let rejection = (deviation > options.fit_threshold).then_some(FitRejection::FitThreshold);
    let mut curves = vec![ra];
    curves.extend(aa.iter().cloned());
    Ok(Cn2Fit {
        candidate: CandidateFit {
            model_copy_number: 2,
            estimated_copy_number: 2.0,
            absolute_deviation: finite(deviation),
            selected: false,
            rejection,
            curves,
        },
        aa,
    })
}

fn fit_cn3(
    distribution: &ChromosomeDistribution,
    options: PolysomyOptions,
    aa: Option<&FitCurve>,
) -> Result<CandidateFit> {
    let x = x_values(distribution);
    let y = y_values(distribution);
    let xrr = x[distribution.rr_boundary];
    let xra = x[distribution.heterozygous_center];
    let xaa = x[distribution.aa_boundary];
    let minimum_separation = 0.5 - 1.0 / (options.minimum_fraction + 2.0);
    let left = Peak::bounded_gaussian(1.0, 1.0 / 3.0, 0.03, xrr, xra - minimum_separation, 7)?
        .scan(1, xrr, xra - minimum_separation)?;
    let right = Peak::bounded_gaussian(1.0, 2.0 / 3.0, 0.03, xra + minimum_separation, xaa, 7)?
        .scan(1, xra + minimum_separation, xaa)?;
    let preliminary = fit(
        vec![left, right],
        &x[distribution.rr_boundary..=distribution.aa_boundary],
        &y[distribution.rr_boundary..=distribution.aa_boundary],
    )?;
    let left = preliminary.peaks[0].physical_parameters();
    let right = preliminary.peaks[1].physical_parameters();
    let separation = ((0.5 - left[1]) + (right[1] - 0.5)) * 0.5;
    let separation = separation.min(0.5 / 3.0);
    let outcome = fit(
        vec![
            Peak::gaussian(left[0], 0.5 - separation, left[2], 5),
            Peak::gaussian(right[0], 0.5 + separation, right[2], 5),
        ],
        &x[distribution.rr_boundary..=distribution.aa_boundary],
        &y[distribution.rr_boundary..=distribution.aa_boundary],
    )?;
    let left = outcome.peaks[0].physical_parameters();
    let right = outcome.peaks[1].physical_parameters();
    let symmetry = ratio(left[0].powi(2), right[0].powi(2));
    let estimated_copy_number = 2.0 + (1.0 - 2.0 * left[1]) / left[1];
    let widths_valid = [left[2], right[2]]
        .into_iter()
        .all(|width| (0.01..=0.3).contains(&width));
    let ra = curve(outcome, distribution.rr_boundary, distribution.aa_boundary);
    let deviation = ra.absolute_deviation + aa.map_or(0.0, |curve| curve.absolute_deviation);
    let rejection = if !widths_valid || deviation > options.fit_threshold {
        Some(FitRejection::FitThreshold)
    } else if symmetry < options.peak_symmetry {
        Some(FitRejection::PeakSymmetry)
    } else {
        None
    };
    let mut curves = vec![ra];
    if let Some(aa) = aa {
        curves.push(aa.clone());
    }
    Ok(CandidateFit {
        model_copy_number: 3,
        estimated_copy_number,
        absolute_deviation: finite(deviation),
        selected: false,
        rejection,
        curves,
    })
}

fn fit_cn4(
    distribution: &ChromosomeDistribution,
    options: PolysomyOptions,
) -> Result<CandidateFit> {
    let x = x_values(distribution);
    let y = y_values(distribution);
    let xrr = x[distribution.rr_boundary];
    let xra = x[distribution.heterozygous_center];
    let xaa = x[distribution.aa_boundary];
    let xmax = x[distribution.fitted_end];
    let aa = if options.include_aa {
        let exponential = Peak::exponential(0.5, 1.0, 0.2, 5).scan(2, 0.01, 0.3)?;
        let gaussian = Peak::bounded_gaussian(0.4, (xaa + xmax) * 0.5, 0.02, xaa, xmax, 7)?
            .scan(1, xaa, xmax)?;
        let outcome = fit(
            vec![exponential, gaussian],
            &x[distribution.aa_boundary..=distribution.fitted_end],
            &y[distribution.aa_boundary..=distribution.fitted_end],
        )?;
        Some(curve(
            outcome,
            distribution.aa_boundary,
            distribution.fitted_end,
        ))
    } else {
        None
    };
    let minimum_separation = 0.25 * options.minimum_fraction;
    let center = Peak::gaussian(1.0, 0.5, 0.03, 5);
    let left = Peak::bounded_gaussian(0.6, 0.3, 0.03, xrr, xra - minimum_separation, 7)?.scan(
        2,
        xrr,
        xra - minimum_separation,
    )?;
    let preliminary = fit(
        vec![center, left],
        &x[distribution.rr_boundary..=distribution.heterozygous_center],
        &y[distribution.rr_boundary..=distribution.heterozygous_center],
    )?;
    let center = preliminary.peaks[0].physical_parameters();
    let left = preliminary.peaks[1].physical_parameters();
    let separation = (0.5 - left[1]).min(0.25);
    let upper = center[0].max(0.100_001);
    let right = Peak::gaussian(left[0], 0.5 + separation, left[2], 5)
        .scan(0, 0.1, upper)?
        .scan(2, 0.01, 0.1)?;
    let outcome = fit(
        vec![
            Peak::gaussian(center[0], 0.5, center[2], 5),
            Peak::gaussian(left[0], 0.5 - separation, left[2], 5),
            right,
        ],
        &x[distribution.rr_boundary..=distribution.aa_boundary],
        &y[distribution.rr_boundary..=distribution.aa_boundary],
    )?;
    let center = outcome.peaks[0].physical_parameters();
    let left = outcome.peaks[1].physical_parameters();
    let right = outcome.peaks[2].physical_parameters();
    let center_size = if center[0] == 0.0 {
        f64::INFINITY
    } else {
        center[0].powi(2)
    };
    let left_size = left[0].powi(2);
    let right_size = right[0].powi(2);
    let left_ratio = ratio(left_size, center_size);
    let right_ratio = ratio(right_size, center_size);
    let symmetry = ratio(left_ratio, right_ratio);
    let minimum_size = left_size.min(right_size) / center_size;
    let placement = (right[1] - 0.5) - (0.5 - left[1]);
    let estimated_copy_number = 3.0 + right[1] - left[1];
    let widths_valid = [center[2], left[2], right[2]]
        .into_iter()
        .all(|width| (0.01..=0.3).contains(&width));
    let ra = curve(outcome, distribution.rr_boundary, distribution.aa_boundary);
    let deviation =
        ra.absolute_deviation + aa.as_ref().map_or(0.0, |curve| curve.absolute_deviation);
    let rejection = if !widths_valid || deviation > options.fit_threshold {
        Some(FitRejection::FitThreshold)
    } else if minimum_size < options.minimum_peak_size {
        Some(FitRejection::PeakSize)
    } else if symmetry < options.peak_symmetry {
        Some(FitRejection::PeakSymmetry)
    } else if placement > 0.1 {
        Some(FitRejection::PeakPlacement)
    } else {
        None
    };
    let mut curves = vec![ra];
    curves.extend(aa);
    Ok(CandidateFit {
        model_copy_number: 4,
        estimated_copy_number,
        absolute_deviation: finite(deviation),
        selected: false,
        rejection,
        curves,
    })
}

fn curve(outcome: FitOutcome, start_bin: usize, end_bin: usize) -> FitCurve {
    FitCurve {
        start_bin,
        end_bin,
        absolute_deviation: outcome.deviation,
        function: outcome
            .peaks
            .iter()
            .map(Peak::function)
            .collect::<Vec<_>>()
            .join(" + "),
    }
}

fn x_values(distribution: &ChromosomeDistribution) -> Vec<f64> {
    distribution.bins.iter().map(|bin| bin.baf).collect()
}

fn y_values(distribution: &ChromosomeDistribution) -> Vec<f64> {
    distribution
        .bins
        .iter()
        .map(|bin| bin.normalized_count)
        .collect()
}

fn finite(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}

fn ratio(first: f64, second: f64) -> f64 {
    if first < second {
        first / second
    } else {
        second / first
    }
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
