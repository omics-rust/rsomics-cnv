use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::allele_frequency::AlleleFrequencyReader;
use crate::emission::{AlleleFrequencies, EmissionModel, EvidenceParameters, SampleParameters};
use crate::hmm::{Hmm, Inference, phred_error_probability};
use crate::selection::{CompiledSelection, SiteSelection};
use crate::signals::{Measurement, RequiredSignals, SampleSelection, SelectedSignalReader};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CallOptions {
    pub sample: SampleParameters,
    pub control_sample: Option<SampleParameters>,
    pub evidence: EvidenceParameters,
    pub transition_probability: f64,
    pub same_state_probability: f64,
    pub lrr_smoothing_window: usize,
    pub optimize_aberrant_fraction: Option<f64>,
}

impl Default for CallOptions {
    fn default() -> Self {
        Self {
            sample: SampleParameters::default(),
            control_sample: None,
            evidence: EvidenceParameters::default(),
            transition_probability: 1e-9,
            same_state_probability: 0.5,
            lrr_smoothing_window: 10,
            optimize_aberrant_fraction: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct CallResult {
    pub sample: String,
    pub control_sample: Option<String>,
    pub chromosomes: Vec<ChromosomeCall>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ChromosomeCall {
    pub reference_name: String,
    pub sites: Vec<SiteCall>,
    pub regions: Vec<RegionCall>,
    pub query_estimate: Option<AberrantEstimate>,
    pub control_estimate: Option<AberrantEstimate>,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct AberrantEstimate {
    pub fraction: f64,
    pub baf_deviation: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SiteCall {
    pub position: u32,
    pub copy_number: u8,
    pub control_copy_number: Option<u8>,
    pub state_probability: f64,
    pub posterior: [f64; 4],
    pub control_posterior: Option<[f64; 4]>,
    pub measurement: Measurement,
    pub modeled_lrr: Option<f64>,
    pub control_measurement: Option<Measurement>,
    pub control_modeled_lrr: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
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
    query: SampleParameters,
    control: Option<SampleParameters>,
    evidence: EvidenceParameters,
    single_hmm: Hmm,
    paired_hmm: Option<Hmm>,
    smoothing_window: usize,
    optimization_minimum: Option<f64>,
}

#[derive(Default)]
struct ChromosomeBuffer {
    positions: Vec<u32>,
    query: Vec<Measurement>,
    control: Vec<Measurement>,
    frequencies: Vec<AlleleFrequencies>,
}

pub fn analyze(
    input: &Path,
    selection: SampleSelection,
    options: CallOptions,
) -> Result<CallResult> {
    let engine = CallEngine::new(options, selection.control.is_some())?;
    analyze_inner(input, selection, engine, None, &SiteSelection::default())
}

/// Runs copy-number inference after applying indexed regions or streaming targets.
pub fn analyze_selected(
    input: &Path,
    selection: SampleSelection,
    options: CallOptions,
    sites: &SiteSelection,
) -> Result<CallResult> {
    let engine = CallEngine::new(options, selection.control.is_some())?;
    analyze_inner(input, selection, engine, None, sites)
}

pub fn analyze_with_allele_frequencies(
    input: &Path,
    selection: SampleSelection,
    options: CallOptions,
    allele_frequencies: &Path,
) -> Result<CallResult> {
    let engine = CallEngine::new(options, selection.control.is_some())?;
    analyze_inner(
        input,
        selection,
        engine,
        Some(allele_frequencies),
        &SiteSelection::default(),
    )
}

/// Runs selected-site inference restricted by an external allele-frequency table.
pub fn analyze_with_allele_frequencies_selected(
    input: &Path,
    selection: SampleSelection,
    options: CallOptions,
    allele_frequencies: &Path,
    sites: &SiteSelection,
) -> Result<CallResult> {
    let engine = CallEngine::new(options, selection.control.is_some())?;
    analyze_inner(input, selection, engine, Some(allele_frequencies), sites)
}

fn analyze_inner(
    input: &Path,
    selection: SampleSelection,
    engine: CallEngine,
    allele_frequency_path: Option<&Path>,
    site_selection: &SiteSelection,
) -> Result<CallResult> {
    let required = if engine.evidence.lrr_weight == 0.0 {
        RequiredSignals::Baf
    } else {
        RequiredSignals::BafAndLrr
    };
    let site_selection = CompiledSelection::new(site_selection)?;
    let mut reader = SelectedSignalReader::open(input, selection, required, site_selection)?;
    let mut allele_frequencies = allele_frequency_path
        .map(|path| AlleleFrequencyReader::open(path, &reader.reference_names()))
        .transpose()?;
    let sample = reader.query_sample().to_owned();
    let control_sample = reader.control_sample().map(str::to_owned);
    let mut chromosomes = Vec::new();
    let mut current_name = None;
    let mut buffer = ChromosomeBuffer::default();

    reader.visit(allele_frequencies.is_some(), |site| {
        let frequencies = if let Some(reader) = &mut allele_frequencies {
            let Some(frequencies) = reader.frequencies(
                site.reference,
                site.position,
                site.alleles.as_deref().unwrap(),
            )?
            else {
                return Ok(());
            };
            frequencies
        } else {
            AlleleFrequencies::default()
        };
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
        buffer.frequencies.push(frequencies);
        if let Some(measurement) = site.control {
            buffer.control.push(measurement);
        }
        Ok(())
    })?;
    if let Some(reference_name) = current_name {
        chromosomes.push(call_chromosome(reference_name, &mut buffer, &engine)?);
    }
    if let Some(reader) = allele_frequencies {
        reader.finish()?;
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

impl CallEngine {
    fn new(options: CallOptions, paired: bool) -> Result<Self> {
        if options.lrr_smoothing_window == 0 {
            return Err(invalid("LRR smoothing window must be positive"));
        }
        if options
            .optimize_aberrant_fraction
            .is_some_and(|minimum| !minimum.is_finite() || !(0.0..=1.0).contains(&minimum))
        {
            return Err(invalid(
                "optimization minimum must be finite and between 0 and 1",
            ));
        }
        let control = match (paired, options.control_sample) {
            (true, value) => Some(value.unwrap_or(options.sample)),
            (false, Some(_)) => {
                return Err(invalid(
                    "control sample parameters require a control sample",
                ));
            }
            (false, None) => None,
        };
        EmissionModel::new(options.sample, options.evidence)?;
        if let Some(parameters) = control {
            EmissionModel::new(parameters, options.evidence)?;
        }
        Ok(Self {
            query: options.sample,
            control,
            evidence: options.evidence,
            single_hmm: Hmm::single_sample(options.transition_probability)?,
            paired_hmm: paired
                .then(|| {
                    Hmm::paired_samples(
                        options.transition_probability,
                        options.same_state_probability,
                    )
                })
                .transpose()?,
            smoothing_window: options.lrr_smoothing_window,
            optimization_minimum: options
                .optimize_aberrant_fraction
                .filter(|minimum| *minimum > 0.0),
        })
    }
}

fn call_chromosome(
    reference_name: String,
    buffer: &mut ChromosomeBuffer,
    engine: &CallEngine,
) -> Result<ChromosomeCall> {
    let query_lrr = smooth_lrr(&buffer.query, engine.smoothing_window);
    let control_lrr = if engine.paired_hmm.is_some() {
        if buffer.control.len() != buffer.query.len() {
            return Err(invalid("query and control site counts differ"));
        }
        Some(smooth_lrr(&buffer.control, engine.smoothing_window))
    } else {
        None
    };

    let (query_parameters, control_parameters) =
        optimize(buffer, &query_lrr, control_lrr.as_deref(), engine)?;
    let inference = infer(
        buffer,
        &query_lrr,
        control_lrr.as_deref(),
        engine,
        query_parameters,
        control_parameters,
    )?;
    let path = inference.path;
    let posterior = inference.posterior;
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
            modeled_lrr: query_lrr[index],
            control_measurement: buffer.control.get(index).copied(),
            control_modeled_lrr: control_lrr.as_ref().and_then(|values| values[index]),
        });
    }
    let regions = regions(&sites)?;
    buffer.positions.clear();
    buffer.query.clear();
    buffer.control.clear();
    buffer.frequencies.clear();
    Ok(ChromosomeCall {
        reference_name,
        sites,
        regions,
        query_estimate: engine
            .optimization_minimum
            .map(|_| AberrantEstimate::from(query_parameters)),
        control_estimate: engine
            .optimization_minimum
            .zip(control_parameters)
            .map(|(_, parameters)| AberrantEstimate::from(parameters)),
    })
}

impl From<SampleParameters> for AberrantEstimate {
    fn from(parameters: SampleParameters) -> Self {
        Self {
            fraction: parameters.aberrant_fraction,
            baf_deviation: parameters.baf_deviation,
        }
    }
}

fn infer(
    buffer: &ChromosomeBuffer,
    query_lrr: &[Option<f64>],
    control_lrr: Option<&[Option<f64>]>,
    engine: &CallEngine,
    query_parameters: SampleParameters,
    control_parameters: Option<SampleParameters>,
) -> Result<Inference> {
    let query_model = EmissionModel::new(query_parameters, engine.evidence)?;
    let query_emissions = emissions(&query_model, &buffer.query, query_lrr, &buffer.frequencies)?;
    if let Some(hmm) = &engine.paired_hmm {
        let control_parameters = control_parameters
            .ok_or_else(|| invalid("paired inference has no control parameters"))?;
        let control_lrr =
            control_lrr.ok_or_else(|| invalid("paired inference has no control LRR values"))?;
        let control_model = EmissionModel::new(control_parameters, engine.evidence)?;
        let control_emissions = emissions(
            &control_model,
            &buffer.control,
            control_lrr,
            &buffer.frequencies,
        )?;
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
        hmm.infer(&buffer.positions, &emissions)
    } else {
        engine.single_hmm.infer(&buffer.positions, &query_emissions)
    }
}

fn emissions(
    model: &EmissionModel,
    measurements: &[Measurement],
    lrr: &[Option<f64>],
    frequencies: &[AlleleFrequencies],
) -> Result<Vec<[f64; 4]>> {
    measurements
        .iter()
        .copied()
        .zip(lrr.iter().copied())
        .zip(frequencies.iter().copied())
        .map(|((measurement, lrr), frequencies)| {
            model.probabilities(Measurement { lrr, ..measurement }, frequencies)
        })
        .collect()
}

fn optimize(
    buffer: &ChromosomeBuffer,
    query_lrr: &[Option<f64>],
    control_lrr: Option<&[Option<f64>]>,
    engine: &CallEngine,
) -> Result<(SampleParameters, Option<SampleParameters>)> {
    let Some(minimum) = engine.optimization_minimum else {
        return Ok((engine.query, engine.control));
    };
    let mut query = engine.query;
    let mut control = engine.control;
    for iteration in 0..20 {
        let inference = infer(buffer, query_lrr, control_lrr, engine, query, control)?;
        let query_converged = update_parameters(
            &buffer.query,
            &inference.posterior,
            false,
            engine.paired_hmm.is_some(),
            engine.query,
            &mut query,
            minimum,
        )?;
        let control_converged = if let Some(parameters) = &mut control {
            update_parameters(
                &buffer.control,
                &inference.posterior,
                true,
                true,
                engine
                    .control
                    .ok_or_else(|| invalid("paired optimization has no control parameters"))?,
                parameters,
                minimum,
            )?
        } else {
            true
        };
        if query_converged && control_converged {
            return Ok((query, control));
        }
        if iteration == 19 {
            return Ok((engine.query, engine.control));
        }
    }
    unreachable!()
}

fn update_parameters(
    measurements: &[Measurement],
    posterior: &[Vec<f64>],
    control: bool,
    paired: bool,
    defaults: SampleParameters,
    parameters: &mut SampleParameters,
    minimum: f64,
) -> Result<bool> {
    if measurements.len() != posterior.len() {
        return Err(invalid(
            "optimization posterior count differs from site count",
        ));
    }
    let mut values = Vec::new();
    let mut aa_variance = 0.0;
    let mut aa_count = 0usize;
    for (measurement, probabilities) in measurements.iter().zip(posterior) {
        let Some(mut baf) = measurement.baf else {
            continue;
        };
        if baf > 0.8 {
            aa_variance += (1.0 - baf).powi(2);
            aa_count += 1;
            continue;
        }
        if baf > 0.5 {
            baf = 1.0 - baf;
        }
        if baf < 0.2 {
            continue;
        }
        values.push((baf, cn3_probability(probabilities, control, paired)?));
    }
    if values.is_empty() {
        parameters.aberrant_fraction = 1.0;
        return Ok(true);
    }
    let smoothed = smooth_values(
        &values
            .iter()
            .map(|(_, probability)| *probability as f32)
            .collect::<Vec<_>>(),
        50,
    );
    let weight = smoothed.iter().map(|value| f64::from(*value)).sum::<f64>();
    if weight == 0.0 {
        parameters.aberrant_fraction = 1.0;
        return Ok(true);
    }
    let mean = values
        .iter()
        .zip(&smoothed)
        .map(|((baf, _), probability)| baf * f64::from(*probability))
        .sum::<f64>()
        / weight;
    let mut variance = values
        .iter()
        .zip(&smoothed)
        .map(|((baf, _), probability)| (baf - mean).powi(2) * f64::from(*probability))
        .sum::<f64>()
        / weight;
    if aa_count > 0 {
        variance = variance.max(aa_variance / aa_count as f64);
    }
    let maximum_mean = 0.5 - variance.sqrt() * 1.644_854;
    if !maximum_mean.is_finite() || maximum_mean <= 0.0 {
        return Err(invalid(
            "aberrant-fraction optimization produced an invalid detection boundary",
        ));
    }
    let mut fraction = mean.recip() - 2.0;
    if mean > maximum_mean || fraction < minimum {
        parameters.aberrant_fraction = 1.0;
        return Ok(true);
    }
    fraction = fraction.min(1.0);
    let converged = (fraction - parameters.aberrant_fraction).abs() < 0.1;
    let default_variance = defaults.baf_deviation.powi(2);
    variance = variance.clamp(0.5 * default_variance, 3.0 * default_variance);
    parameters.aberrant_fraction = fraction;
    parameters.baf_deviation = variance.sqrt();
    Ok(converged)
}

fn cn3_probability(probabilities: &[f64], control: bool, paired: bool) -> Result<f64> {
    if !paired {
        return probabilities.get(3).copied().ok_or_else(|| {
            invalid("single-sample optimization posterior has fewer than 4 states")
        });
    }
    if probabilities.len() != 16 {
        return Err(invalid(format!(
            "paired optimization posterior has {} states instead of 16",
            probabilities.len()
        )));
    }
    Ok(if control {
        (0..4).map(|query| probabilities[query * 4 + 3]).sum()
    } else {
        probabilities[12..16].iter().sum()
    })
}

fn smooth_lrr(measurements: &[Measurement], window: usize) -> Vec<Option<f64>> {
    if window <= 1 || measurements.is_empty() {
        return measurements
            .iter()
            .map(|measurement| measurement.lrr)
            .collect();
    }
    let values = measurements
        .iter()
        .map(|measurement| measurement.lrr.unwrap_or(0.0) as f32)
        .collect::<Vec<_>>();
    let smoothed = smooth_values(&values, window);
    measurements
        .iter()
        .zip(smoothed)
        .map(|(measurement, value)| {
            if measurement.baf.is_some() {
                Some(f64::from(value))
            } else {
                measurement.lrr
            }
        })
        .collect()
}

fn smooth_values(values: &[f32], window: usize) -> Vec<f32> {
    if window <= 1 || values.is_empty() {
        return values.to_vec();
    }
    let mut smoothed = vec![0.0f32; values.len()];
    let left = window / 2;
    let right = window - left;
    for (index, output) in smoothed.iter_mut().enumerate() {
        let start = index.saturating_sub(left);
        let end = (index + right).min(values.len());
        let sum = values[start..end].iter().copied().sum::<f32>();
        *output = sum / (end - start) as f32;
    }
    smoothed
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
