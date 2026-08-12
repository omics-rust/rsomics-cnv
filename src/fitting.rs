use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{DMatrix, DVector, Dyn, storage::Owned};
use rsomics_common::{Result, RsomicsError};

const PARAMETERS: usize = 3;

#[derive(Clone, Copy)]
enum PeakKind {
    Gaussian,
    BoundedGaussian { minimum: f64, maximum: f64 },
    Exponential,
}

#[derive(Clone, Copy)]
struct Scan {
    minimum: f64,
    maximum: f64,
    iterations: usize,
}

#[derive(Clone)]
pub(crate) struct Peak {
    kind: PeakKind,
    parameters: [f64; PARAMETERS],
    original: [f64; PARAMETERS],
    fitted: [bool; PARAMETERS],
    scans: [Option<Scan>; PARAMETERS],
}

impl Peak {
    pub(crate) fn gaussian(scale: f64, center: f64, deviation: f64, mask: u8) -> Self {
        Self::new(PeakKind::Gaussian, scale, center, deviation, mask)
    }

    pub(crate) fn bounded_gaussian(
        scale: f64,
        center: f64,
        deviation: f64,
        minimum: f64,
        maximum: f64,
        mask: u8,
    ) -> Result<Self> {
        if !minimum.is_finite() || !maximum.is_finite() || minimum >= maximum {
            return Err(invalid(
                "bounded Gaussian requires an increasing finite interval",
            ));
        }
        let mut peak = Self::new(
            PeakKind::BoundedGaussian { minimum, maximum },
            scale,
            center,
            deviation,
            mask,
        );
        peak.set_physical(1, center);
        peak.original = peak.parameters;
        Ok(peak)
    }

    pub(crate) fn exponential(scale: f64, center: f64, deviation: f64, mask: u8) -> Self {
        Self::new(PeakKind::Exponential, scale, center, deviation, mask)
    }

    fn new(kind: PeakKind, scale: f64, center: f64, deviation: f64, mask: u8) -> Self {
        Self {
            kind,
            parameters: [scale, center, deviation],
            original: [scale, center, deviation],
            fitted: std::array::from_fn(|index| mask & (1 << index) != 0),
            scans: [None; PARAMETERS],
        }
    }

    pub(crate) fn scan(mut self, parameter: usize, minimum: f64, maximum: f64) -> Result<Self> {
        if parameter >= PARAMETERS
            || !minimum.is_finite()
            || !maximum.is_finite()
            || minimum >= maximum
        {
            return Err(invalid("invalid peak scan interval"));
        }
        self.scans[parameter] = Some(Scan {
            minimum,
            maximum,
            iterations: 16,
        });
        Ok(self)
    }

    pub(crate) fn physical_parameters(&self) -> [f64; PARAMETERS] {
        let mut parameters = self.parameters.map(f64::abs);
        if let PeakKind::BoundedGaussian { minimum, maximum } = self.kind {
            parameters[1] = 0.5 * (self.parameters[1].cos() + 1.0) * (maximum - minimum) + minimum;
        }
        parameters
    }

    pub(crate) fn function(&self) -> String {
        let [scale, center, deviation] = self.physical_parameters();
        match self.kind {
            PeakKind::Exponential => {
                format!("{scale:.6}**2 * exp((x-{center:.6})/{deviation:.6}**2)")
            }
            _ => format!("{scale:.6}**2 * exp(-(x-{center:.6})**2/{deviation:.6}**2)"),
        }
    }

    fn set_physical(&mut self, parameter: usize, value: f64) {
        self.parameters[parameter] = match (self.kind, parameter) {
            (PeakKind::BoundedGaussian { minimum, maximum }, 1) => {
                let value = value.clamp(minimum, maximum);
                (2.0 * (value - minimum) / (maximum - minimum) - 1.0).acos()
            }
            _ => value,
        };
    }

    pub(crate) fn value(&self, x: f64) -> f64 {
        let [scale, center, deviation] = self.physical_parameters();
        match self.kind {
            PeakKind::Exponential => scale * scale * ((x - center) / deviation.powi(2)).exp(),
            _ => scale * scale * (-((x - center) / deviation).powi(2)).exp(),
        }
    }

    fn derivative(&self, parameter: usize, x: f64) -> f64 {
        let scale = self.parameters[0];
        let deviation = self.parameters[2];
        match self.kind {
            PeakKind::Gaussian => {
                let center = self.parameters[1];
                gaussian_derivative(scale, center, deviation, parameter, x)
            }
            PeakKind::BoundedGaussian { minimum, maximum } => {
                let angle = self.parameters[1];
                let center = 0.5 * (angle.cos() + 1.0) * (maximum - minimum) + minimum;
                let offset = x - center;
                let exponential = (-(offset / deviation).powi(2)).exp();
                match parameter {
                    0 => 2.0 * scale * exponential,
                    1 => {
                        -scale.powi(2) * angle.sin() * (maximum - minimum) * offset * exponential
                            / deviation.powi(2)
                    }
                    2 => 2.0 * scale.powi(2) * offset.powi(2) * exponential / deviation.powi(3),
                    _ => 0.0,
                }
            }
            PeakKind::Exponential => {
                let center = self.parameters[1];
                let exponential = ((x - center) / deviation.powi(2)).exp();
                match parameter {
                    0 => 2.0 * scale * exponential,
                    2 => -2.0 * scale.powi(2) * (x - center) * exponential / deviation.powi(3),
                    _ => 0.0,
                }
            }
        }
    }
}

fn gaussian_derivative(scale: f64, center: f64, deviation: f64, parameter: usize, x: f64) -> f64 {
    let offset = x - center;
    let exponential = (-(offset / deviation).powi(2)).exp();
    match parameter {
        0 => 2.0 * scale * exponential,
        1 => 2.0 * scale.powi(2) * offset * exponential / deviation.powi(2),
        2 => 2.0 * scale.powi(2) * offset.powi(2) * exponential / deviation.powi(3),
        _ => 0.0,
    }
}

struct FitProblem {
    peaks: Vec<Peak>,
    variables: Vec<(usize, usize)>,
    parameters: DVector<f64>,
    x: Vec<f64>,
    y: Vec<f64>,
}

impl FitProblem {
    fn new(peaks: Vec<Peak>, x: &[f64], y: &[f64]) -> Self {
        let variables = peaks
            .iter()
            .enumerate()
            .flat_map(|(peak, value)| {
                value
                    .fitted
                    .iter()
                    .enumerate()
                    .filter_map(move |(parameter, fitted)| fitted.then_some((peak, parameter)))
            })
            .collect::<Vec<_>>();
        let parameters = DVector::from_iterator(
            variables.len(),
            variables
                .iter()
                .map(|&(peak, parameter)| peaks[peak].parameters[parameter]),
        );
        Self {
            peaks,
            variables,
            parameters,
            x: x.to_vec(),
            y: y.to_vec(),
        }
    }

    fn finite(&self) -> bool {
        self.peaks
            .iter()
            .flat_map(|peak| peak.parameters)
            .all(f64::is_finite)
    }

    fn deviation(&self) -> f64 {
        self.x
            .iter()
            .zip(&self.y)
            .map(|(&x, &y)| (self.peaks.iter().map(|peak| peak.value(x)).sum::<f64>() - y).abs())
            .sum()
    }
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for FitProblem {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, parameters: &DVector<f64>) {
        self.parameters.copy_from(parameters);
        for (value, &(peak, parameter)) in parameters.iter().zip(&self.variables) {
            self.peaks[peak].parameters[parameter] = *value;
        }
    }

    fn params(&self) -> DVector<f64> {
        self.parameters.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let residuals = DVector::from_iterator(
            self.x.len(),
            self.x.iter().zip(&self.y).map(|(&x, &y)| {
                (self.peaks.iter().map(|peak| peak.value(x)).sum::<f64>() - y) / 0.01
            }),
        );
        residuals
            .iter()
            .all(|value| value.is_finite())
            .then_some(residuals)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let mut jacobian = DMatrix::zeros(self.x.len(), self.variables.len());
        for (column, &(peak, parameter)) in self.variables.iter().enumerate() {
            for (row, &x) in self.x.iter().enumerate() {
                jacobian[(row, column)] = self.peaks[peak].derivative(parameter, x) / 0.01;
            }
        }
        jacobian
            .iter()
            .all(|value| value.is_finite())
            .then_some(jacobian)
    }
}

pub(crate) struct FitOutcome {
    pub(crate) peaks: Vec<Peak>,
    pub(crate) deviation: f64,
}

pub(crate) fn fit(peaks: Vec<Peak>, x: &[f64], y: &[f64]) -> Result<FitOutcome> {
    if x.len() != y.len() || x.is_empty() {
        return Err(invalid(
            "peak fitting requires equal nonempty x and y vectors",
        ));
    }
    let variables = peaks
        .iter()
        .map(|peak| peak.fitted.iter().filter(|fitted| **fitted).count())
        .sum::<usize>();
    if variables == 0 || x.len() < variables {
        return Err(invalid(
            "peak fitting has an invalid observation-to-parameter ratio",
        ));
    }
    let iterations = peaks
        .iter()
        .flat_map(|peak| peak.scans.iter().flatten())
        .map(|scan| scan.iterations)
        .max()
        .unwrap_or(0);
    let mut random = Random::new(0x6a09_e667_f3bc_c909);
    let mut best = None;
    for _ in 0..=iterations {
        let mut trial = peaks.clone();
        for peak in &mut trial {
            peak.parameters = peak.original;
            for parameter in 0..PARAMETERS {
                if let Some(scan) = peak.scans[parameter] {
                    let value = scan.minimum + random.next() * (scan.maximum - scan.minimum);
                    peak.set_physical(parameter, value);
                }
            }
        }
        let problem = FitProblem::new(trial, x, y);
        let (problem, _) = LevenbergMarquardt::new()
            .with_tol(1e-8)
            .with_patience(150)
            .minimize(problem);
        if !problem.finite() {
            continue;
        }
        let deviation = problem.deviation();
        if !deviation.is_finite() {
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|best: &FitOutcome| deviation < best.deviation)
        {
            best = Some(FitOutcome {
                peaks: problem.peaks,
                deviation,
            });
        }
    }
    best.ok_or_else(|| invalid("peak fitting produced no finite solution"))
}

struct Random(u64);

impl Random {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let value = self.0.wrapping_mul(0x2545_f491_4f6c_dd1d);
        (value >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analytic_derivatives_match_finite_differences() {
        let peaks = [
            Peak::gaussian(0.8, 0.47, 0.06, 7),
            Peak::bounded_gaussian(0.8, 0.47, 0.06, 0.3, 0.7, 7).unwrap(),
            Peak::exponential(0.8, 1.0, 0.2, 5),
        ];
        for peak in peaks {
            for parameter in 0..PARAMETERS {
                if !peak.fitted[parameter] {
                    continue;
                }
                let mut lower = peak.clone();
                let mut upper = peak.clone();
                lower.parameters[parameter] -= 1e-7;
                upper.parameters[parameter] += 1e-7;
                let numerical = (upper.value(0.43) - lower.value(0.43)) / 2e-7;
                let analytic = peak.derivative(parameter, 0.43);
                assert!(
                    (numerical - analytic).abs() <= 1e-6 * analytic.abs().max(1.0),
                    "{numerical} != {analytic}"
                );
            }
        }
    }

    #[test]
    fn full_jacobian_recovers_a_gaussian() {
        let expected = Peak::gaussian(0.9, 0.51, 0.045, 7);
        let x = (0..150)
            .map(|index| index as f64 / 149.0)
            .collect::<Vec<_>>();
        let y = x
            .iter()
            .map(|&value| expected.value(value))
            .collect::<Vec<_>>();
        let outcome = fit(vec![Peak::gaussian(0.7, 0.45, 0.08, 7)], &x, &y).unwrap();
        let actual = outcome.peaks[0].physical_parameters();
        for (actual, expected) in actual.into_iter().zip([0.9, 0.51, 0.045]) {
            assert!((actual - expected).abs() < 1e-6, "{actual} != {expected}");
        }
        assert!(outcome.deviation < 2e-6, "{}", outcome.deviation);
    }
}
