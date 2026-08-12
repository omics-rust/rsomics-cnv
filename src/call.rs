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
    pub same_state_probability: f64,
    pub lrr_smoothing_window: usize,
}

impl Default for CallOptions {
    fn default() -> Self {
        Self {
            sample: SampleParameters::default(),
            evidence: EvidenceParameters::default(),
            transition_probability: 1e-9,
            same_state_probability: 0.5,
            lrr_smoothing_window: 10,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CallResult {
    pub sample: String,
    pub control_sample: Option<String>,
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
    pub control_copy_number: Option<u8>,
    pub state_probability: f64,
    pub posterior: [f64; 4],
    pub control_posterior: Option<[f64; 4]>,
    pub measurement: Measurement,
    pub control_measurement: Option<Measurement>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RegionCall {
    pub start: u32,
    pub end: u32,
    pub copy_number: u8,
    pub control_copy_number: Option<u8>,
    pub quality: f64,
    pub sites: usize,
    pub heterozygous_sites: usize,
    pub control_heterozygous_sites: Option<usize>,
}

struct CallEngine {
    model: EmissionModel,
    single_hmm: Hmm,
    paired_hmm: Option<Hmm>,
    smoothing_window: usize,
}

#[derive(Default)]
struct ChromosomeBuffer {
    positions: Vec<u32>,
    query: Vec<Measurement>,
    control: Vec<Measurement>,
}

pub fn analyze(
    input: &Path,
    selection: SampleSelection,
    options: CallOptions,
) -> Result<CallResult> {
    if options.lrr_smoothing_window == 0 {
        return Err(invalid("LRR smoothing window must be positive"));
    }
    let engine = CallEngine {
        model: EmissionModel::new(options.sample, options.evidence)?,
        single_hmm: Hmm::single_sample(options.transition_probability)?,
        paired_hmm: selection
            .control
            .as_ref()
            .map(|_| {
                Hmm::paired_samples(
                    options.transition_probability,
                    options.same_state_probability,
                )
            })
            .transpose()?,
        smoothing_window: options.lrr_smoothing_window,
    };
    let required = if options.evidence.lrr_weight == 0.0 {
        RequiredSignals::Baf
    } else {
        RequiredSignals::BafAndLrr
    };
    let mut reader = SignalReader::open(input, selection, required)?;
    let sample = reader.query_sample().to_owned();
    let control_sample = reader.control_sample().map(str::to_owned);
    let mut chromosomes = Vec::new();
    let mut current_name = None;
    let mut buffer = ChromosomeBuffer::default();

    while let Some(site) = reader.next_site()? {
        if current_name
            .as_deref()
            .is_some_and(|name| name != site.reference_name)
        {
            chromosomes.push(call_chromosome(
                current_name.take().unwrap(),
                &mut buffer,
                &engine,
            )?);
        }
        current_name.get_or_insert(site.reference_name);
        buffer.positions.push(site.position);
        buffer.query.push(site.query);
        if let Some(measurement) = site.control {
            buffer.control.push(measurement);
        }
    }
    if let Some(reference_name) = current_name {
        chromosomes.push(call_chromosome(reference_name, &mut buffer, &engine)?);
    }
    if chromosomes.is_empty() {
        return Err(invalid(
            "no informative BAF sites remain after sample selection",
        ));
    }

    Ok(CallResult {
        sample,
        control_sample,
        chromosomes,
    })
}

fn call_chromosome(
    reference_name: String,
    buffer: &mut ChromosomeBuffer,
    engine: &CallEngine,
) -> Result<ChromosomeCall> {
    smooth_lrr(&mut buffer.query, engine.smoothing_window)?;
    if engine.paired_hmm.is_some() {
        if buffer.control.len() != buffer.query.len() {
            return Err(invalid("query and control site counts differ"));
        }
        smooth_lrr(&mut buffer.control, engine.smoothing_window)?;
    }

    let query_emissions = buffer
        .query
        .iter()
        .copied()
        .map(|measurement| {
            engine
                .model
                .probabilities(measurement, AlleleFrequencies::default())
        })
        .collect::<Result<Vec<_>>>()?;
    let (path, posterior) = if let Some(hmm) = &engine.paired_hmm {
        let control_emissions = buffer
            .control
            .iter()
            .copied()
            .map(|measurement| {
                engine
                    .model
                    .probabilities(measurement, AlleleFrequencies::default())
            })
            .collect::<Result<Vec<_>>>()?;
        let emissions = query_emissions
            .iter()
            .zip(control_emissions)
            .map(|(query, control)| {
                let mut joint = [0.0; 16];
                for query_state in 0..4 {
                    for control_state in 0..4 {
                        joint[query_state * 4 + control_state] =
                            query[query_state] * control[control_state];
                    }
                }
                joint
            })
            .collect::<Vec<_>>();
        let inference = hmm.infer(&buffer.positions, &emissions)?;
        (inference.path, inference.posterior)
    } else {
        let inference = engine
            .single_hmm
            .infer(&buffer.positions, &query_emissions)?;
        (inference.path, inference.posterior)
    };
    let mut sites = Vec::with_capacity(buffer.positions.len());
    for (index, ((position, measurement), joint_posterior)) in buffer
        .positions
        .iter()
        .copied()
        .zip(buffer.query.iter().copied())
        .zip(posterior)
        .enumerate()
    {
        let (copy_number, control_copy_number, state_probability, posterior, control_posterior) =
            if engine.paired_hmm.is_some() {
                let state = path[index];
                let mut query = [0.0; 4];
                let mut control = [0.0; 4];
                for query_state in 0..4 {
                    for control_state in 0..4 {
                        let probability = joint_posterior[query_state * 4 + control_state];
                        query[query_state] += probability;
                        control[control_state] += probability;
                    }
                }
                (
                    u8::try_from(state / 4)
                        .map_err(|_| invalid("query copy-number state exceeds uint8"))?,
                    Some(
                        u8::try_from(state % 4)
                            .map_err(|_| invalid("control copy-number state exceeds uint8"))?,
                    ),
                    joint_posterior[state],
                    query,
                    Some(control),
                )
            } else {
                let posterior: [f64; 4] =
                    joint_posterior.try_into().map_err(|values: Vec<f64>| {
                        invalid(format!(
                            "HMM returned {} posterior states instead of 4",
                            values.len()
                        ))
                    })?;
                (
                    u8::try_from(path[index])
                        .map_err(|_| invalid("copy-number state exceeds uint8"))?,
                    None,
                    posterior[path[index]],
                    posterior,
                    None,
                )
            };
        sites.push(SiteCall {
            position,
            copy_number,
            control_copy_number,
            state_probability,
            posterior,
            control_posterior,
            measurement,
            control_measurement: buffer.control.get(index).copied(),
        });
    }
    let regions = regions(&sites)?;
    buffer.positions.clear();
    buffer.query.clear();
    buffer.control.clear();
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
        .map(|measurement| measurement.lrr.unwrap_or(0.0) as f32)
        .collect::<Vec<_>>();
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
        if measurement.baf.is_some() {
            measurement.lrr = Some(f64::from(value));
        }
    }
    Ok(())
}

fn regions(sites: &[SiteCall]) -> Result<Vec<RegionCall>> {
    let mut regions = Vec::new();
    let mut start = 0;
    for index in 1..=sites.len() {
        if index < sites.len()
            && sites[index].copy_number == sites[start].copy_number
            && sites[index].control_copy_number == sites[start].control_copy_number
        {
            continue;
        }
        let slice = &sites[start..index];
        let copy_number = sites[start].copy_number;
        let control_copy_number = sites[start].control_copy_number;
        let mean =
            slice.iter().map(|site| site.state_probability).sum::<f64>() / slice.len() as f64;
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
            control_copy_number,
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
            control_heterozygous_sites: control_copy_number.map(|_| {
                slice
                    .iter()
                    .filter(|site| {
                        site.control_measurement
                            .and_then(|measurement| measurement.baf)
                            .is_some_and(|baf| baf > 0.25 && baf < 0.75)
                    })
                    .count()
            }),
        });
        start = index;
    }
    Ok(regions)
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
