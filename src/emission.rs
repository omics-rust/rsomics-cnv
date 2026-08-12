use std::f64::consts::{PI, SQRT_2};

use rsomics_common::{Result, RsomicsError};

use crate::signals::Measurement;

const CN0: usize = 0;
const CN1: usize = 1;
const CN2: usize = 2;
const CN3: usize = 3;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SampleParameters {
    pub baf_deviation: f64,
    pub lrr_deviation: f64,
    pub aberrant_fraction: f64,
}

impl Default for SampleParameters {
    fn default() -> Self {
        Self {
            baf_deviation: 0.04,
            lrr_deviation: 0.2,
            aberrant_fraction: 1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvidenceParameters {
    pub baf_weight: f64,
    pub lrr_weight: f64,
    pub error_probability: f64,
}

impl Default for EvidenceParameters {
    fn default() -> Self {
        Self {
            baf_weight: 1.0,
            lrr_weight: 0.2,
            error_probability: 1e-4,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlleleFrequencies {
    pub rr: f64,
    pub ra: f64,
    pub aa: f64,
}

impl Default for AlleleFrequencies {
    fn default() -> Self {
        Self {
            rr: 0.76,
            ra: 0.14,
            aa: 0.098,
        }
    }
}

impl AlleleFrequencies {
    pub fn new(rr: f64, ra: f64, aa: f64) -> Result<Self> {
        let values = [rr, ra, aa];
        if values
            .iter()
            .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
        {
            return Err(invalid(
                "genotype frequencies must be finite and between 0 and 1",
            ));
        }
        let sum: f64 = values.iter().sum();
        if sum <= 0.0 {
            return Err(invalid("genotype frequencies have zero mass"));
        }
        Ok(Self { rr, ra, aa })
    }

    pub fn from_non_reference_frequency(frequency: f64) -> Result<Self> {
        if !frequency.is_finite() || !(0.0..=1.0).contains(&frequency) {
            return Err(invalid(
                "non-reference allele frequency must be finite and between 0 and 1",
            ));
        }
        let reference = 1.0 - frequency;
        Self::new(
            reference * reference,
            2.0 * reference * frequency,
            frequency * frequency,
        )
    }
}

#[derive(Clone, Copy, Debug)]
struct Peak {
    mean: f64,
    normalization: f64,
}

#[derive(Clone, Debug)]
pub struct EmissionModel {
    sample: SampleParameters,
    evidence: EvidenceParameters,
    peaks: [Peak; 9],
}

impl EmissionModel {
    pub fn new(sample: SampleParameters, evidence: EvidenceParameters) -> Result<Self> {
        validate_parameters(sample, evidence)?;
        let fraction = sample.aberrant_fraction;
        let means = [
            0.0,
            1.0,
            0.0,
            0.5,
            1.0,
            0.0,
            1.0 / (2.0 + fraction),
            (1.0 + fraction) / (2.0 + fraction),
            1.0,
        ];
        let peaks = means.map(|mean| Peak {
            mean,
            normalization: normal_mass(mean, sample.baf_deviation),
        });
        if peaks
            .iter()
            .any(|peak| !peak.normalization.is_finite() || peak.normalization <= 0.0)
        {
            return Err(invalid(
                "BAF deviation produces an empty or non-finite Gaussian peak",
            ));
        }
        Ok(Self {
            sample,
            evidence,
            peaks,
        })
    }

    pub fn probabilities(
        &self,
        measurement: Measurement,
        frequencies: AlleleFrequencies,
    ) -> Result<[f64; 4]> {
        let frequencies = AlleleFrequencies::new(frequencies.rr, frequencies.ra, frequencies.aa)?;
        let Some(baf) = measurement.baf else {
            if measurement.lrr.is_some() {
                return Err(invalid("LRR is present while BAF is missing"));
            }
            return Ok([0.5, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]);
        };
        if !baf.is_finite() || !(0.0..=1.0).contains(&baf) {
            return Err(invalid("BAF must be finite and between 0 and 1"));
        }
        let lrr = match measurement.lrr {
            Some(value) if value.is_finite() => value,
            Some(_) => return Err(invalid("LRR must be finite")),
            None if self.evidence.lrr_weight == 0.0 => 0.0,
            None => {
                return Err(invalid(
                    "LRR is required when its evidence weight is nonzero",
                ));
            }
        };

        let density = |index| self.density(baf, index);
        let cn1 = density(0) * (frequencies.rr + frequencies.ra * 0.5)
            + density(1) * (frequencies.aa + frequencies.ra * 0.5);
        let cn2 =
            density(2) * frequencies.rr + density(3) * frequencies.ra + density(4) * frequencies.aa;
        let cn3 = density(5) * frequencies.rr
            + density(6) * frequencies.ra * 0.5
            + density(7) * frequencies.ra * 0.5
            + density(8) * frequencies.aa;
        let total = cn1 + cn2 + cn3;
        if !total.is_finite() || total <= 0.0 {
            return Err(invalid(
                "BAF peaks have zero or non-finite probability mass",
            ));
        }
        let baf = [cn1 / total, cn2 / total, cn3 / total];

        let lrr_variance = self.sample.lrr_deviation * self.sample.lrr_deviation;
        let lrr = [
            (-(lrr + 0.45).powi(2) / lrr_variance).exp(),
            (-(lrr - 0.00).powi(2) / lrr_variance).exp(),
            (-(lrr - 0.30).powi(2) / lrr_variance).exp(),
        ];
        let mut output = [0.0; 4];
        output[CN0] = 0.0;
        for (index, state) in [CN1, CN2, CN3].into_iter().enumerate() {
            output[state] = self.evidence.error_probability
                + (1.0 - self.evidence.baf_weight + self.evidence.baf_weight * baf[index])
                    * (1.0 - self.evidence.lrr_weight + self.evidence.lrr_weight * lrr[index]);
        }
        if output
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
            || output.iter().all(|value| *value == 0.0)
        {
            return Err(invalid(
                "emission calculation produced invalid probability mass",
            ));
        }
        Ok(output)
    }

    pub(crate) fn lrr_weight(&self) -> f64 {
        self.evidence.lrr_weight
    }

    fn density(&self, baf: f64, peak: usize) -> f64 {
        let peak = self.peaks[peak];
        let variance = self.sample.baf_deviation * self.sample.baf_deviation;
        (-(baf - peak.mean).powi(2) * 0.5 / variance).exp()
            / peak.normalization
            / (2.0 * PI * variance).sqrt()
    }
}

fn validate_parameters(sample: SampleParameters, evidence: EvidenceParameters) -> Result<()> {
    for (name, value) in [
        ("BAF deviation", sample.baf_deviation),
        ("LRR deviation", sample.lrr_deviation),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid(format!("{name} must be finite and positive")));
        }
    }
    for (name, value) in [
        ("aberrant fraction", sample.aberrant_fraction),
        ("BAF weight", evidence.baf_weight),
        ("LRR weight", evidence.lrr_weight),
        ("error probability", evidence.error_probability),
    ] {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(invalid(format!(
                "{name} must be finite and between 0 and 1"
            )));
        }
    }
    Ok(())
}

fn normal_mass(mean: f64, deviation: f64) -> f64 {
    let top = 1.0 - 0.5 * libm::erfc((1.0 - mean) / (deviation * SQRT_2));
    let bottom = 1.0 - 0.5 * libm::erfc(-mean / (deviation * SQRT_2));
    top - bottom
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
