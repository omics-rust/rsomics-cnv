use std::path::{Path, PathBuf};
use std::process;
use std::str::FromStr;

use clap::{Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{OutputArgs, Result, RsomicsError, ToolMeta, run as run_tool};
use serde::Serialize;

use crate::call::{
    CallOptions, analyze_selected as analyze_calls, analyze_with_allele_frequencies_selected,
};
use crate::emission::{EvidenceParameters, SampleParameters};
use crate::polysomy::{PolysomyOptions, analyze_selected as analyze_polysomy};
use crate::reports::{
    CallReportOptions, PolysomyReportOptions, write_call_reports_with_options,
    write_polysomy_reports_with_options,
};
use crate::selection::{OverlapMode, SiteSelection};
use crate::signals::SampleSelection;

const META: ToolMeta = ToolMeta {
    name: "rsomics-cnv",
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-cnv",
    version,
    about = "BAF and LRR copy-number analysis workflows",
    arg_required_else_help = true,
    subcommand_required = true
)]
struct Cli {
    #[command(flatten)]
    output: OutputArgs,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Infer site and region copy-number states with an HMM
    Call(CallArgs),
    /// Estimate chromosome copy number from BAF peak distributions
    Polysomy(PolysomyArgs),
}

#[derive(Debug, Args)]
struct CallArgs {
    /// Coordinate-sorted VCF or BCF containing FORMAT/BAF and FORMAT/LRR
    #[arg(value_name = "VCF_OR_BCF")]
    input: PathBuf,

    /// New directory for compatibility and machine-readable reports
    #[arg(short, long, value_name = "DIRECTORY")]
    output: PathBuf,

    /// Query sample name; defaults to the only sample in the input
    #[arg(short, long, value_name = "NAME")]
    sample: Option<String>,

    /// Matched control sample name
    #[arg(short, long, value_name = "NAME")]
    control: Option<String>,

    /// Tab-separated CHROM, POS, REF,ALT, AF file; restricts analysis to listed sites
    #[arg(short = 'f', long = "allele-frequencies", value_name = "TSV")]
    allele_frequencies: Option<PathBuf>,

    #[command(flatten)]
    selection: SelectionArgs,

    /// Query BAF Gaussian deviation; append ,CONTROL for a matched control
    #[arg(
        short = 'd',
        long,
        default_value = "0.04",
        value_name = "QUERY[,CONTROL]"
    )]
    baf_deviation: SampleParameter,

    /// Query LRR Gaussian deviation; append ,CONTROL for a matched control
    #[arg(
        short = 'k',
        long,
        default_value = "0.2",
        value_name = "QUERY[,CONTROL]"
    )]
    lrr_deviation: SampleParameter,

    /// Query aberrant-cell fraction; append ,CONTROL for a matched control
    #[arg(
        short = 'a',
        long,
        default_value = "1.0",
        value_name = "QUERY[,CONTROL]"
    )]
    aberrant_fraction: SampleParameter,

    /// Estimate aberrant-cell fraction per chromosome down to this minimum
    #[arg(short = 'O', long = "optimize", value_name = "FRACTION")]
    optimize_aberrant_fraction: Option<f64>,

    /// Relative contribution of BAF evidence
    #[arg(short = 'b', long, default_value_t = 1.0)]
    baf_weight: f64,

    /// Relative contribution of LRR evidence; zero allows BAF-only input
    #[arg(short = 'l', long, default_value_t = 0.2)]
    lrr_weight: f64,

    /// Emission-model error floor
    #[arg(short = 'e', long, default_value_t = 1e-4)]
    error_probability: f64,

    /// HMM transition probability per base
    #[arg(short = 'x', long, default_value_t = 1e-9)]
    transition_probability: f64,

    /// Prior probability that query and control share a state
    #[arg(short = 'P', long, default_value_t = 0.5)]
    same_state_probability: f64,

    /// Number of neighboring sites used for LRR smoothing
    #[arg(short = 'L', long, default_value_t = 10)]
    lrr_smoothing_window: usize,

    /// Write chromosome SVG plots whose maximum region quality reaches this value
    #[arg(short = 'p', long, value_name = "QUALITY")]
    plot_threshold: Option<f64>,
}

#[derive(Debug, Args)]
struct PolysomyArgs {
    /// Coordinate-sorted VCF or BCF containing FORMAT/BAF
    #[arg(value_name = "VCF_OR_BCF")]
    input: PathBuf,

    /// New directory for compatibility and machine-readable reports
    #[arg(short, long, value_name = "DIRECTORY")]
    output: PathBuf,

    /// Sample name; defaults to the only sample in the input
    #[arg(short, long, value_name = "NAME")]
    sample: Option<String>,

    #[command(flatten)]
    selection: SelectionArgs,

    /// Maximum accepted absolute fit deviation
    #[arg(short = 'f', long, default_value_t = 3.3)]
    fit_threshold: f64,

    /// Improvement required before selecting a higher copy-number model
    #[arg(short = 'c', long, default_value_t = 0.7)]
    copy_number_penalty: f64,

    /// Minimum symmetry of paired BAF peaks
    #[arg(short = 'p', long, default_value_t = 0.5)]
    peak_symmetry: f64,

    /// Minimum peak height relative to the fitted distribution
    #[arg(short = 'b', long, default_value_t = 0.1)]
    minimum_peak_size: f64,

    /// Minimum chromosome fraction used for a preliminary CN1 decision
    #[arg(short = 'm', long, default_value_t = 0.1)]
    minimum_fraction: f64,

    /// Include the homozygous alternate peak in model fitting
    #[arg(short = 'i', long)]
    include_aa: bool,

    /// Distribution smoothing parameter; negative values use increasing windows
    #[arg(long, default_value_t = -3, allow_hyphen_values = true)]
    smoothing: i32,

    /// Write BAF-distribution and chromosome copy-number SVG plots
    #[arg(long)]
    plots: bool,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum OverlapArg {
    #[value(name = "0")]
    Position,
    #[value(name = "1")]
    Record,
    #[value(name = "2")]
    Variant,
}

#[derive(Clone, Copy, Debug)]
struct SampleParameter {
    query: f64,
    control: Option<f64>,
}

impl FromStr for SampleParameter {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        let mut values = value.split(',');
        let query = parse_parameter(values.next().unwrap_or_default())?;
        let control = values.next().map(parse_parameter).transpose()?;
        if values.next().is_some() {
            return Err("expected QUERY or QUERY,CONTROL".to_owned());
        }
        Ok(Self { query, control })
    }
}

fn parse_parameter(value: &str) -> std::result::Result<f64, String> {
    value
        .parse()
        .map_err(|_| format!("{value:?} is not a floating-point value"))
}

impl From<OverlapArg> for OverlapMode {
    fn from(value: OverlapArg) -> Self {
        match value {
            OverlapArg::Position => Self::Position,
            OverlapArg::Record => Self::Record,
            OverlapArg::Variant => Self::Variant,
        }
    }
}

#[derive(Debug, Args)]
struct SelectionArgs {
    /// Comma-separated indexed regions
    #[arg(
        short = 'r',
        long,
        value_delimiter = ',',
        conflicts_with = "regions_file"
    )]
    regions: Vec<String>,

    /// BED, VCF, or tabular indexed regions
    #[arg(short = 'R', long, value_name = "FILE", conflicts_with = "regions")]
    regions_file: Option<PathBuf>,

    /// Region inclusion rule: POS, record span, or variant span
    #[arg(long, default_value = "1")]
    regions_overlap: OverlapArg,

    /// Comma-separated streaming targets; prefix the first target with ^ to exclude
    #[arg(
        short = 't',
        long,
        value_delimiter = ',',
        conflicts_with = "targets_file"
    )]
    targets: Vec<String>,

    /// BED, VCF, or tabular streaming targets; prefix the path with ^ to exclude
    #[arg(short = 'T', long, value_name = "FILE", conflicts_with = "targets")]
    targets_file: Option<PathBuf>,

    /// Target inclusion rule: POS, record span, or variant span
    #[arg(long, default_value = "0")]
    targets_overlap: OverlapArg,
}

impl SelectionArgs {
    fn build(mut self) -> Result<SiteSelection> {
        let mut exclude_targets = false;
        if let Some(first) = self.targets.first_mut()
            && let Some(target) = first.strip_prefix('^')
        {
            if target.is_empty() {
                return Err(RsomicsError::ConfigError(
                    "target list is empty after ^".to_owned(),
                ));
            }
            *first = target.to_owned();
            exclude_targets = true;
        }
        if let Some(path) = &self.targets_file
            && let Some(value) = path.to_str().and_then(|value| value.strip_prefix('^'))
        {
            if value.is_empty() {
                return Err(RsomicsError::ConfigError(
                    "target file path is empty after ^".to_owned(),
                ));
            }
            self.targets_file = Some(PathBuf::from(value));
            exclude_targets = true;
        }
        Ok(SiteSelection {
            regions: self.regions,
            regions_file: self.regions_file,
            regions_overlap: self.regions_overlap.into(),
            targets: self.targets,
            targets_file: self.targets_file,
            targets_overlap: self.targets_overlap.into(),
            exclude_targets,
        })
    }
}

#[derive(Debug, Serialize)]
struct RunSummary {
    workflow: &'static str,
    sample: String,
    output: PathBuf,
    chromosomes: usize,
    sites: Option<usize>,
    regions: Option<usize>,
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli.command))
}

fn execute(command: Command) -> Result<RunSummary> {
    match command {
        Command::Call(args) => {
            require_new_output(&args.output)?;
            let (sample, control_sample) = call_sample_parameters(&args)?;
            let site_selection = args.selection.build()?;
            let selection = SampleSelection {
                query: args.sample,
                control: args.control,
            };
            let options = CallOptions {
                sample,
                control_sample,
                evidence: EvidenceParameters {
                    baf_weight: args.baf_weight,
                    lrr_weight: args.lrr_weight,
                    error_probability: args.error_probability,
                },
                transition_probability: args.transition_probability,
                same_state_probability: args.same_state_probability,
                lrr_smoothing_window: args.lrr_smoothing_window,
                optimize_aberrant_fraction: args.optimize_aberrant_fraction,
            };
            let result = if let Some(frequencies) = args.allele_frequencies {
                analyze_with_allele_frequencies_selected(
                    &args.input,
                    selection,
                    options,
                    &frequencies,
                    &site_selection,
                )?
            } else {
                analyze_calls(&args.input, selection, options, &site_selection)?
            };
            write_call_reports_with_options(
                &args.output,
                &result,
                CallReportOptions {
                    plot_threshold: args.plot_threshold,
                },
            )?;
            Ok(RunSummary {
                workflow: "call",
                sample: result.sample,
                output: args.output,
                chromosomes: result.chromosomes.len(),
                sites: Some(
                    result
                        .chromosomes
                        .iter()
                        .map(|chromosome| chromosome.sites.len())
                        .sum(),
                ),
                regions: Some(
                    result
                        .chromosomes
                        .iter()
                        .map(|chromosome| chromosome.regions.len())
                        .sum(),
                ),
            })
        }
        Command::Polysomy(args) => {
            require_new_output(&args.output)?;
            let site_selection = args.selection.build()?;
            let result = analyze_polysomy(
                &args.input,
                args.sample,
                PolysomyOptions {
                    fit_threshold: args.fit_threshold,
                    copy_number_penalty: args.copy_number_penalty,
                    peak_symmetry: args.peak_symmetry,
                    minimum_peak_size: args.minimum_peak_size,
                    minimum_fraction: args.minimum_fraction,
                    include_aa: args.include_aa,
                    smoothing: args.smoothing,
                },
                &site_selection,
            )?;
            write_polysomy_reports_with_options(
                &args.output,
                &result,
                PolysomyReportOptions { plots: args.plots },
            )?;
            Ok(RunSummary {
                workflow: "polysomy",
                sample: result.sample,
                output: args.output,
                chromosomes: result.chromosomes.len(),
                sites: None,
                regions: None,
            })
        }
    }
}

fn call_sample_parameters(args: &CallArgs) -> Result<(SampleParameters, Option<SampleParameters>)> {
    let sample = SampleParameters {
        baf_deviation: args.baf_deviation.query,
        lrr_deviation: args.lrr_deviation.query,
        aberrant_fraction: args.aberrant_fraction.query,
    };
    let distinct = [
        args.baf_deviation.control,
        args.lrr_deviation.control,
        args.aberrant_fraction.control,
    ]
    .into_iter()
    .any(|value| value.is_some());
    if args.control.is_none() && distinct {
        return Err(RsomicsError::InvalidInput(
            "paired model parameters require --control".to_owned(),
        ));
    }
    let control = args.control.as_ref().map(|_| SampleParameters {
        baf_deviation: args
            .baf_deviation
            .control
            .unwrap_or(args.baf_deviation.query),
        lrr_deviation: args
            .lrr_deviation
            .control
            .unwrap_or(args.lrr_deviation.query),
        aberrant_fraction: args
            .aberrant_fraction
            .control
            .unwrap_or(args.aberrant_fraction.query),
    });
    Ok((sample, control))
}

fn require_new_output(output: &Path) -> Result<()> {
    if output.exists() {
        return Err(RsomicsError::ConfigError(format!(
            "output directory {} already exists",
            output.display()
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    #[test]
    fn command_tree_is_valid() {
        rsomics_help::command::<Cli>().debug_assert();
    }

    #[test]
    fn commands_share_the_product_help_tree() {
        let top = Cli::command().render_long_help().to_string();
        assert!(top.contains("call"), "{top}");
        assert!(top.contains("polysomy"), "{top}");
        let error = Cli::try_parse_from(["rsomics-cnv", "polysomy", "--help"]).unwrap_err();
        let help = error.to_string();
        assert!(help.contains("--fit-threshold"), "{help}");
        assert!(help.contains("--include-aa"), "{help}");
    }
}
