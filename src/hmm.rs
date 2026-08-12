use std::collections::HashMap;

use rsomics_common::{Result, RsomicsError};

const COPY_NUMBER_STATES: usize = 4;
const MAX_EXACT_POWER: usize = 10_000;

#[derive(Clone, Debug)]
pub struct Hmm {
    states: usize,
    transition: Vec<f64>,
    initial: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Inference {
    pub path: Vec<usize>,
    pub posterior: Vec<Vec<f64>>,
}

impl Hmm {
    pub fn single_sample(jump_probability: f64) -> Result<Self> {
        validate_jump(jump_probability)?;
        let states = COPY_NUMBER_STATES;
        let stay = 1.0 - jump_probability * (states - 1) as f64;
        let mut transition = vec![jump_probability; states * states];
        for state in 0..states {
            transition[state * states + state] = stay;
        }
        let initial = vec![1.0 / 6.0, 1.0 / 6.0, 0.5, 1.0 / 6.0];
        Ok(Self {
            states,
            transition,
            initial,
        })
    }

    pub fn paired_samples(jump_probability: f64, same_probability: f64) -> Result<Self> {
        validate_jump(jump_probability)?;
        if !same_probability.is_finite() || !(0.0..=1.0).contains(&same_probability) {
            return Err(invalid(
                "same-state probability must be finite and between 0 and 1",
            ));
        }

        let states = COPY_NUMBER_STATES * COPY_NUMBER_STATES;
        let stay = 1.0 - jump_probability * (COPY_NUMBER_STATES - 1) as f64;
        let jump = (1.0 - stay) / (states - 1) as f64;
        let mut transition = vec![0.0; states * states];
        for source in 0..states {
            let (source_query, source_control) = paired_state(source);
            let mut sum = 0.0;
            for destination in 0..states {
                let (query, control) = paired_state(destination);
                let query_probability = if query == source_query { stay } else { jump };
                let control_probability = if control == source_control {
                    stay
                } else {
                    jump
                };
                let independent = query_probability * control_probability;
                let probability = if query == control && source_query == source_control {
                    independent * (1.0 - same_probability) + independent.sqrt() * same_probability
                } else if query == control {
                    independent
                } else {
                    independent * (1.0 - same_probability)
                };
                transition[destination * states + source] = probability;
                sum += probability;
            }
            if !sum.is_finite() || sum <= 0.0 {
                return Err(invalid(
                    "paired transition parameters produce an empty probability column",
                ));
            }
            for destination in 0..states {
                transition[destination * states + source] /= sum;
            }
        }

        let single = [1.0 / 6.0, 1.0 / 6.0, 0.5, 1.0 / 6.0];
        let mut initial = vec![0.0; states];
        for (state, probability) in initial.iter_mut().enumerate() {
            let (query, control) = paired_state(state);
            *probability = single[query] * single[control];
            if query != control {
                *probability *= 1.0 - same_probability;
            }
        }
        normalize(&mut initial, "paired initial probabilities")?;

        Ok(Self {
            states,
            transition,
            initial,
        })
    }

    pub fn infer<const N: usize>(
        &self,
        positions: &[u32],
        emissions: &[[f64; N]],
    ) -> Result<Inference> {
        if N != self.states {
            return Err(invalid(format!(
                "expected {} emission probabilities per site, found {N}",
                self.states
            )));
        }
        if positions.len() != emissions.len() {
            return Err(invalid(format!(
                "position count {} does not match emission count {}",
                positions.len(),
                emissions.len()
            )));
        }
        if positions.is_empty() {
            return Err(invalid("HMM input has no sites"));
        }
        if positions.windows(2).any(|pair| pair[1] < pair[0]) {
            return Err(invalid("HMM positions must be sorted"));
        }
        for (site, values) in emissions.iter().enumerate() {
            if values
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
            {
                return Err(invalid(format!(
                    "site {} has a non-finite or negative emission probability",
                    site + 1
                )));
            }
            if values.iter().all(|value| *value == 0.0) {
                return Err(invalid(format!(
                    "site {} has no positive emission probability",
                    site + 1
                )));
            }
        }

        let mut transitions =
            TransitionCache::new(self.states, &self.transition, maximum_step(positions));
        let (path, forward) = self.forward_viterbi(positions, emissions, &mut transitions)?;
        let posterior = self.backward(positions, emissions, forward, &mut transitions)?;
        Ok(Inference { path, posterior })
    }

    fn forward_viterbi<const N: usize>(
        &self,
        positions: &[u32],
        emissions: &[[f64; N]],
        transitions: &mut TransitionCache,
    ) -> Result<(Vec<usize>, Vec<Vec<f64>>)> {
        let count = positions.len();
        let mut trace = vec![vec![0usize; self.states]; count];
        let mut viterbi = self.initial.clone();
        let mut forward = Vec::with_capacity(count + 1);
        forward.push(self.initial.clone());

        for site in 0..count {
            let steps = step_count(positions, site);
            let transition = transitions.for_steps(steps);
            let mut next_viterbi = vec![0.0; self.states];
            let mut next_forward = vec![0.0; self.states];
            for destination in 0..self.states {
                let row = destination * self.states;
                let mut best_source = 0;
                let mut best = viterbi[0] * transition[row];
                let mut total = forward[site][0] * transition[row];
                for source in 1..self.states {
                    let probability = viterbi[source] * transition[row + source];
                    if probability > best {
                        best = probability;
                        best_source = source;
                    }
                    total += forward[site][source] * transition[row + source];
                }
                trace[site][destination] = best_source;
                next_viterbi[destination] = best * emissions[site][destination];
                next_forward[destination] = total * emissions[site][destination];
            }
            normalize(
                &mut next_viterbi,
                &format!("Viterbi probabilities at site {}", site + 1),
            )?;
            normalize(
                &mut next_forward,
                &format!("forward probabilities at site {}", site + 1),
            )?;
            viterbi = next_viterbi;
            forward.push(next_forward);
        }

        let mut final_state = 0;
        for state in 1..self.states {
            if viterbi[state] > viterbi[final_state] {
                final_state = state;
            }
        }
        let mut path = vec![0; count];
        let mut state = final_state;
        for site in (0..count).rev() {
            path[site] = state;
            state = trace[site][state];
        }
        Ok((path, forward))
    }

    fn backward<const N: usize>(
        &self,
        positions: &[u32],
        emissions: &[[f64; N]],
        forward: Vec<Vec<f64>>,
        transitions: &mut TransitionCache,
    ) -> Result<Vec<Vec<f64>>> {
        let count = positions.len();
        let mut posterior = vec![vec![0.0; self.states]; count];
        let mut backward = vec![1.0; self.states];
        let mut previous_position = positions[count - 1];

        for site in (0..count).rev() {
            for state in 0..self.states {
                posterior[site][state] = forward[site + 1][state] * backward[state];
            }
            normalize(
                &mut posterior[site],
                &format!("posterior probabilities at site {}", site + 1),
            )?;

            let steps = if positions[site] == previous_position {
                1
            } else {
                previous_position - positions[site]
            };
            previous_position = positions[site];
            let transition = transitions.for_steps(steps);
            let mut next = vec![0.0; self.states];
            for source in 0..self.states {
                for destination in 0..self.states {
                    next[source] += backward[destination]
                        * emissions[site][destination]
                        * transition[destination * self.states + source];
                }
            }
            normalize(
                &mut next,
                &format!("backward probabilities at site {}", site + 1),
            )?;
            backward = next;
        }
        Ok(posterior)
    }
}

pub fn phred_error_probability(probability: f64) -> Result<f64> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(invalid(
            "error probability must be finite and between 0 and 1",
        ));
    }
    if probability == 0.0 {
        Ok(99.0)
    } else {
        Ok((-4.3429 * probability.ln()).min(99.0))
    }
}

struct TransitionCache {
    states: usize,
    powers: Vec<Vec<f64>>,
    large: HashMap<u32, Vec<f64>>,
}

impl TransitionCache {
    fn new(states: usize, base: &[f64], maximum_step: u32) -> Self {
        let count = usize::try_from(maximum_step)
            .unwrap_or(usize::MAX)
            .clamp(1, MAX_EXACT_POWER);
        let mut powers = Vec::with_capacity(count);
        powers.push(base.to_vec());
        while powers.len() < count {
            powers.push(multiply(base, powers.last().unwrap(), states));
        }
        Self {
            states,
            powers,
            large: HashMap::new(),
        }
    }

    fn for_steps(&mut self, steps: u32) -> &[f64] {
        let index = usize::try_from(steps - 1).unwrap_or(usize::MAX);
        if index < self.powers.len() {
            return &self.powers[index];
        }
        if !self.large.contains_key(&steps) {
            let block_count = (steps - 1) / MAX_EXACT_POWER as u32;
            let remainder = usize::try_from((steps - 1) % MAX_EXACT_POWER as u32).unwrap();
            let blocks = matrix_power(&self.powers[MAX_EXACT_POWER - 1], block_count, self.states);
            let value = multiply(&blocks, &self.powers[remainder], self.states);
            self.large.insert(steps, value);
        }
        self.large.get(&steps).unwrap()
    }
}

fn validate_jump(probability: f64) -> Result<()> {
    if !probability.is_finite() || !(0.0..=0.25).contains(&probability) {
        return Err(invalid(
            "transition probability must be finite and between 0 and 0.25",
        ));
    }
    Ok(())
}

fn paired_state(state: usize) -> (usize, usize) {
    (state / COPY_NUMBER_STATES, state % COPY_NUMBER_STATES)
}

fn step_count(positions: &[u32], site: usize) -> u32 {
    if site == 0 || positions[site] == positions[site - 1] {
        1
    } else {
        positions[site] - positions[site - 1]
    }
}

fn maximum_step(positions: &[u32]) -> u32 {
    positions
        .windows(2)
        .map(|pair| pair[1].saturating_sub(pair[0]).max(1))
        .max()
        .unwrap_or(1)
}

fn normalize(values: &mut [f64], name: &str) -> Result<()> {
    let sum: f64 = values.iter().sum();
    if !sum.is_finite() || sum <= 0.0 {
        return Err(invalid(format!("{name} have zero or non-finite mass")));
    }
    for value in values {
        *value /= sum;
    }
    Ok(())
}

fn multiply(left: &[f64], right: &[f64], states: usize) -> Vec<f64> {
    let mut output = vec![0.0; states * states];
    for row in 0..states {
        for column in 0..states {
            for inner in 0..states {
                output[row * states + column] +=
                    left[row * states + inner] * right[inner * states + column];
            }
        }
    }
    output
}

fn matrix_power(matrix: &[f64], mut exponent: u32, states: usize) -> Vec<f64> {
    let mut result = vec![0.0; states * states];
    for state in 0..states {
        result[state * states + state] = 1.0;
    }
    let mut factor = matrix.to_vec();
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = multiply(&factor, &result, states);
        }
        exponent >>= 1;
        if exponent > 0 {
            factor = multiply(&factor, &factor, states);
        }
    }
    result
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
