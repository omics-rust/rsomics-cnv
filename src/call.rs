use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::emission::{AlleleFrequencies, EmissionModel, EvidenceParameters, SampleParameters};
use crate::hmm::{Hmm, phred_error_probability};
use crate::signals::{Measurement, RequiredSignals, SampleSelection, SignalReader};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CallOptions {
    pub sample: SampleParameters,
    pub evidence: EvidenceParameters,
    pub transition_probability: f64,
    pub lrr_smoothing_window: usize,
}

impl Default for CallOptions {
    fn default() -> Self {
        Self {
            sample: SampleParameters::default(),
            evidence: EvidenceParameters::default(),
            transition_probability: 1e-9,
            lrr_smoothing_window: 10,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CallResult {
    pub sample: String,
    pub chromosomes: Vec<ChromosomeCall>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChromosomeCall {
    pub reference_name: String,
    pub sites: Vec<SiteCall>,
    pub regions: Vec<RegionCall>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SiteCall {
    pub position: u32,
    pub copy_number: u8,
    pub posterior: [f64; 4],
    pub measurement: Measurement,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RegionCall {
    pub start: u32,
    pub end: u32,
    pub copy_number: u8,
    pub quality: f64,
    pub sites: usize,
    pub heterozygous_sites: usize,
}

pub fn analyze(
    input: &Path,
    selection: SampleSelection,
    options: CallOptions,
) -> Result<CallResult> {
    if options.lrr_smoothing_window == 0 {
        return Err(invalid("LRR smoothing window must be positive"));
    }
    if selection.control.is_some() {
        return Err(invalid(
            "paired query and control inference is not yet exposed by this analysis entry point",
        ));
    }
    let model = EmissionModel::new(options.sample, options.evidence)?;
    let hmm = Hmm::single_sample(options.transition_probability)?;
    let required = if options.evidence.lrr_weight == 0.0 {
        RequiredSignals::Baf
    } else {
        RequiredSignals::BafAndLrr
    };
    let mut reader = SignalReader::open(input, selection, required)?;
    let sample = reader.query_sample().to_owned();
    let mut chromosomes = Vec::new();
    let mut current_name = None;
    let mut positions = Vec::new();
    let mut measurements = Vec::new();

    while let Some(site) = reader.next_site()? {
        if current_name
            .as_deref()
            .is_some_and(|name| name != site.reference_name)
        {
            chromosomes.push(call_chromosome(
                current_name.take().unwrap(),
                &mut positions,
                &mut measurements,
                &model,
                &hmm,
                options.lrr_smoothing_window,
            )?);
        }
        current_name.get_or_insert(site.reference_name);
        positions.push(site.position);
        measurements.push(site.query);
    }
    if let Some(reference_name) = current_name {
        chromosomes.push(call_chromosome(
            reference_name,
            &mut positions,
            &mut measurements,
            &model,
            &hmm,
            options.lrr_smoothing_window,
        )?);
    }
    if chromosomes.is_empty() {
        return Err(invalid(
            "no informative BAF sites remain after sample selection",
        ));
    }

    Ok(CallResult {
        sample,
        chromosomes,
    })
}

fn call_chromosome(
    reference_name: String,
    positions: &mut Vec<u32>,
    measurements: &mut Vec<Measurement>,
    model: &EmissionModel,
    hmm: &Hmm,
    smoothing_window: usize,
) -> Result<ChromosomeCall> {
    smooth_lrr(measurements, smoothing_window)?;
    let mut emissions = Vec::with_capacity(measurements.len());
    for measurement in measurements.iter().copied() {
        emissions.push(model.probabilities(measurement, AlleleFrequencies::default())?);
    }
    let inference = hmm.infer(positions, &emissions)?;
    let mut sites = Vec::with_capacity(positions.len());
    for (index, ((position, measurement), posterior)) in positions
        .iter()
        .copied()
        .zip(measurements.iter().copied())
        .zip(inference.posterior)
        .enumerate()
    {
        let posterior: [f64; 4] = posterior.try_into().map_err(|values: Vec<f64>| {
            invalid(format!(
                "HMM returned {} posterior states instead of 4",
                values.len()
            ))
        })?;
        sites.push(SiteCall {
            position,
            copy_number: u8::try_from(inference.path[index])
                .map_err(|_| invalid("copy-number state exceeds uint8"))?,
            posterior,
            measurement,
        });
    }
    let regions = regions(&sites)?;
    positions.clear();
    measurements.clear();
    Ok(ChromosomeCall {
        reference_name,
        sites,
        regions,
    })
}

fn smooth_lrr(measurements: &mut [Measurement], window: usize) -> Result<()> {
    if window <= 1 || measurements.is_empty() {
        return Ok(());
    }
    let values = measurements
        .iter()
        .enumerate()
        .map(|(index, measurement)| {
            measurement.lrr.map(|value| value as f32).ok_or_else(|| {
                invalid(format!(
                    "site {} is missing LRR before smoothing",
                    index + 1
                ))
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut smoothed = vec![0.0f32; values.len()];
    let left = window / 2;
    let right = window - left;
    for (index, output) in smoothed.iter_mut().enumerate() {
        let start = index.saturating_sub(left);
        let end = (index + right).min(values.len());
        let sum = values[start..end].iter().copied().sum::<f32>();
        *output = sum / (end - start) as f32;
    }
    for (measurement, value) in measurements.iter_mut().zip(smoothed) {
        measurement.lrr = Some(f64::from(value));
    }
    Ok(())
}

fn regions(sites: &[SiteCall]) -> Result<Vec<RegionCall>> {
    let mut regions = Vec::new();
    let mut start = 0;
    for index in 1..=sites.len() {
        if index < sites.len() && sites[index].copy_number == sites[start].copy_number {
            continue;
        }
        let slice = &sites[start..index];
        let copy_number = sites[start].copy_number;
        let mean = slice
            .iter()
            .map(|site| site.posterior[copy_number as usize])
            .sum::<f64>()
            / slice.len() as f64;
        let error = (1.0 - mean).clamp(0.0, 1.0);
        let end = if index < sites.len() {
            sites[index]
                .position
                .checked_sub(1)
                .ok_or_else(|| invalid("region boundary underflows uint32"))?
        } else {
            sites[index - 1].position
        };
        regions.push(RegionCall {
            start: sites[start].position,
            end,
            copy_number,
            quality: phred_error_probability(error)?,
            sites: slice.len(),
            heterozygous_sites: slice
                .iter()
                .filter(|site| {
                    site.measurement
                        .baf
                        .is_some_and(|baf| baf > 0.25 && baf < 0.75)
                })
                .count(),
        });
        start = index;
    }
    Ok(regions)
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
