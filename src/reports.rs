use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError};
use serde::Serialize;

use crate::call::{CallResult, ChromosomeCall};
use crate::plots::{for_each_call_plot, for_each_polysomy_plot};
use crate::polysomy::{BINS, PolysomyResult};

const DAT_HEADER: &str = "# [1]Chromosome\t[2]Position\t[3]BAF\t[4]LRR\n";
const CN_HEADER: &str =
    "# [1]Chromosome\t[2]Position\t[3]CN\t[4]P(CN0)\t[5]P(CN1)\t[6]P(CN2)\t[7]P(CN3)\n";
const SUMMARY_HEADER: &str = "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy Number state\t[6]Quality\t[7]nSites\t[8]nHETs\n";
const ESTIMATE_HEADER: &str = "# CF, cell fraction estimate\t[2]Chromosome\t[3]Start\t[4]End\t[5]Cell fraction\t[6]BAF deviation\n";

#[derive(Serialize)]
struct ResultDocument<'a, T> {
    schema: &'static str,
    producer: Producer,
    result: &'a T,
}

#[derive(Serialize)]
struct Producer {
    name: &'static str,
    version: &'static str,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CallReportOptions {
    pub plot_threshold: Option<f64>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PolysomyReportOptions {
    pub plots: bool,
}

pub fn write_call_reports(output: &Path, result: &CallResult) -> Result<()> {
    write_call_reports_with_options(output, result, CallReportOptions::default())
}

pub fn write_call_reports_with_options(
    output: &Path,
    result: &CallResult,
    options: CallReportOptions,
) -> Result<()> {
    validate_call_result(result)?;
    write_report_directory(output, |directory| {
        write_sample_reports(directory, result, false)?;
        if result.control_sample.is_some() {
            write_sample_reports(directory, result, true)?;
            write_file(directory.join("summary.tab"), |writer| {
                write_joint_summary(writer, result)
            })?;
        }
        write_json(
            directory.join("result.json"),
            "rsomics-cnv/call-result/v1",
            result,
        )?;
        if let Some(threshold) = options.plot_threshold {
            for_each_call_plot(result, threshold, |name, svg| {
                write_plot(directory, name, svg)
            })?;
        }
        Ok(())
    })
}

pub fn write_polysomy_reports(output: &Path, result: &PolysomyResult) -> Result<()> {
    write_polysomy_reports_with_options(output, result, PolysomyReportOptions::default())
}

pub fn write_polysomy_reports_with_options(
    output: &Path,
    result: &PolysomyResult,
    options: PolysomyReportOptions,
) -> Result<()> {
    validate_polysomy_result(result)?;
    write_report_directory(output, |directory| {
        write_file(directory.join("dist.dat"), |writer| {
            write_polysomy_compatibility(writer, result)
        })?;
        write_json(
            directory.join("result.json"),
            "rsomics-cnv/polysomy-result/v1",
            result,
        )?;
        if options.plots {
            for_each_polysomy_plot(result, |name, svg| write_plot(directory, name, svg))?;
        }
        Ok(())
    })
}

fn write_plot(directory: &Path, name: String, svg: String) -> Result<()> {
    write_file(directory.join(name), |writer| {
        writer.write_all(svg.as_bytes()).map_err(RsomicsError::Io)
    })
}

fn write_report_directory(output: &Path, body: impl FnOnce(&Path) -> Result<()>) -> Result<()> {
    if output.exists() {
        return Err(RsomicsError::ConfigError(format!(
            "output directory {} already exists",
            output.display()
        )));
    }
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let stage = tempfile::Builder::new()
        .prefix(".rsomics-cnv-")
        .tempdir_in(parent)
        .rs_with_context(|| format!("staging report directory beside {}", output.display()))?;
    body(stage.path())?;
    File::open(stage.path())
        .and_then(|directory| directory.sync_all())
        .rs_with_context(|| format!("syncing staged reports for {}", output.display()))?;
    fs::rename(stage.path(), output)
        .rs_with_context(|| format!("committing report directory {}", output.display()))?;
    let _ = stage.keep();
    File::open(parent)
        .and_then(|directory| directory.sync_all())
        .rs_with_context(|| format!("syncing report parent {}", parent.display()))?;
    Ok(())
}

fn write_json<T: Serialize>(path: PathBuf, schema: &'static str, result: &T) -> Result<()> {
    write_file(path, |writer| {
        let document = ResultDocument {
            schema,
            producer: Producer {
                name: "rsomics-cnv",
                version: env!("CARGO_PKG_VERSION"),
            },
            result,
        };
        serde_json::to_writer_pretty(&mut *writer, &document)
            .map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
        writer.write_all(b"\n").map_err(RsomicsError::Io)
    })
}

fn write_polysomy_compatibility(writer: &mut dyn Write, result: &PolysomyResult) -> Result<()> {
    writeln!(
        writer,
        "# This file was produced by: rsomics-cnv {}",
        env!("CARGO_PKG_VERSION")
    )
    .map_err(RsomicsError::Io)?;
    writer.write_all(b"#\n# DIST\t[2]Chrom\t[3]BAF\t[4]Normalized Count\n# FIT\t[2]Chrom\t[3]Goodness of Fit\t[4]iFrom\t[5]iTo\t[6]The Fitted Function\n# CN\t[2]Chrom\t[3]Estimated Copy Number\t[4]Absolute fit deviation\n").map_err(RsomicsError::Io)?;
    for chromosome in &result.chromosomes {
        for bin in &chromosome.distribution.bins {
            writeln!(
                writer,
                "DIST\t{}\t{:.6}\t{:.6}",
                chromosome.distribution.reference_name, bin.baf, bin.normalized_count
            )
            .map_err(RsomicsError::Io)?;
        }
        if let Some(candidate) = chromosome
            .candidates
            .iter()
            .find(|candidate| candidate.selected)
        {
            for curve in &candidate.curves {
                writeln!(
                    writer,
                    "FIT\t{}\t{}\t{}\t{}\t{}",
                    chromosome.distribution.reference_name,
                    scientific(curve.absolute_deviation),
                    curve.start_bin,
                    curve.end_bin,
                    curve.function
                )
                .map_err(RsomicsError::Io)?;
            }
        }
        write!(
            writer,
            "CN\t{}\t{:.2}",
            chromosome.distribution.reference_name, chromosome.copy_number
        )
        .map_err(RsomicsError::Io)?;
        if let Some(deviation) = chromosome.absolute_deviation {
            write!(writer, "\t{deviation:.6}").map_err(RsomicsError::Io)?;
        }
        writeln!(writer).map_err(RsomicsError::Io)?;
    }
    Ok(())
}

fn write_sample_reports(directory: &Path, result: &CallResult, control: bool) -> Result<()> {
    let name = if control {
        result.control_sample.as_deref().ok_or_else(|| {
            RsomicsError::InvalidInput(
                "control report requested without a control sample".to_owned(),
            )
        })?
    } else {
        &result.sample
    };
    write_file(directory.join(format!("dat.{name}.tab")), |writer| {
        write_data(writer, &result.chromosomes, control)
    })?;
    write_file(directory.join(format!("cn.{name}.tab")), |writer| {
        write_copy_number(writer, &result.chromosomes, control)
    })?;
    write_file(directory.join(format!("summary.{name}.tab")), |writer| {
        write_summary(writer, &result.chromosomes, control)
    })
}

fn write_data(writer: &mut dyn Write, chromosomes: &[ChromosomeCall], control: bool) -> Result<()> {
    writer
        .write_all(DAT_HEADER.as_bytes())
        .map_err(RsomicsError::Io)?;
    for chromosome in chromosomes {
        for site in &chromosome.sites {
            let measurement = if control {
                site.control_measurement
                    .ok_or_else(|| inconsistent("missing control measurement"))?
            } else {
                site.measurement
            };
            if let Some(baf) = measurement.baf {
                writeln!(
                    writer,
                    "{}\t{}\t{baf:.3}\t{:.3}",
                    chromosome.reference_name,
                    site.position,
                    measurement.lrr.unwrap_or(0.0)
                )
                .map_err(RsomicsError::Io)?;
            }
        }
    }
    Ok(())
}

fn write_copy_number(
    writer: &mut dyn Write,
    chromosomes: &[ChromosomeCall],
    control: bool,
) -> Result<()> {
    writer
        .write_all(CN_HEADER.as_bytes())
        .map_err(RsomicsError::Io)?;
    for chromosome in chromosomes {
        for site in &chromosome.sites {
            let (copy_number, posterior) = if control {
                (
                    site.control_copy_number
                        .ok_or_else(|| inconsistent("missing control copy-number state"))?,
                    site.control_posterior
                        .ok_or_else(|| inconsistent("missing control posterior"))?,
                )
            } else {
                (site.copy_number, site.posterior)
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                chromosome.reference_name,
                site.position,
                copy_number,
                posterior[0],
                posterior[1],
                posterior[2],
                posterior[3]
            )
            .map_err(RsomicsError::Io)?;
        }
    }
    Ok(())
}

fn write_summary(
    writer: &mut dyn Write,
    chromosomes: &[ChromosomeCall],
    control: bool,
) -> Result<()> {
    writer
        .write_all(SUMMARY_HEADER.as_bytes())
        .map_err(RsomicsError::Io)?;
    if chromosomes.iter().any(|chromosome| {
        if control {
            chromosome.control_estimate.is_some()
        } else {
            chromosome.query_estimate.is_some()
        }
    }) {
        writer
            .write_all(ESTIMATE_HEADER.as_bytes())
            .map_err(RsomicsError::Io)?;
    }
    for chromosome in chromosomes {
        let estimate = if control {
            chromosome.control_estimate
        } else {
            chromosome.query_estimate
        };
        if let Some(estimate) = estimate {
            write_estimate(writer, chromosome, estimate)?;
        }
        for region in &chromosome.regions {
            let copy_number = if control {
                region
                    .control_copy_number
                    .ok_or_else(|| inconsistent("missing control region state"))?
            } else {
                region.copy_number
            };
            let heterozygous_sites = if control {
                region
                    .control_heterozygous_sites
                    .ok_or_else(|| inconsistent("missing control heterozygous-site count"))?
            } else {
                region.heterozygous_sites
            };
            writeln!(
                writer,
                "RG\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
                chromosome.reference_name,
                region.start,
                region.end,
                copy_number,
                region.quality,
                region.sites,
                heterozygous_sites
            )
            .map_err(RsomicsError::Io)?;
        }
    }
    Ok(())
}

fn write_joint_summary(writer: &mut dyn Write, result: &CallResult) -> Result<()> {
    let control = result
        .control_sample
        .as_deref()
        .ok_or_else(|| inconsistent("missing control sample"))?;
    writeln!(
        writer,
        "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:{}\t[6]Copy number:{}\t[7]Quality\t[8]nSites in (5)\t[9]nHETs in (5)\t[10]nSites in (6)\t[11]nHETs in (6)",
        result.sample, control
    )
    .map_err(RsomicsError::Io)?;
    if result
        .chromosomes
        .iter()
        .any(|chromosome| chromosome.query_estimate.is_some())
    {
        writeln!(
            writer,
            "# CF, cell fraction estimate\t[2]Chromosome\t[3]Start\t[4]End\t[5]Cell fraction:{}\t[6]Cell fraction:{}\t[7]BAF deviation:{}\t[8]BAF deviation:{}",
            result.sample, control, result.sample, control
        )
        .map_err(RsomicsError::Io)?;
    }
    for chromosome in &result.chromosomes {
        if let Some(query) = chromosome.query_estimate {
            let control = chromosome
                .control_estimate
                .ok_or_else(|| inconsistent("missing control aberrant-fraction estimate"))?;
            let (start, end) = chromosome_bounds(chromosome)?;
            writeln!(
                writer,
                "CF\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.6}\t{:.6}",
                chromosome.reference_name,
                start,
                end,
                query.fraction,
                control.fraction,
                query.baf_deviation,
                control.baf_deviation
            )
            .map_err(RsomicsError::Io)?;
        }
        for region in &chromosome.regions {
            writeln!(
                writer,
                "RG\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
                chromosome.reference_name,
                region.start,
                region.end,
                region.copy_number,
                region
                    .control_copy_number
                    .ok_or_else(|| inconsistent("missing control region state"))?,
                region.quality,
                region.sites,
                region.heterozygous_sites,
                region.sites,
                region
                    .control_heterozygous_sites
                    .ok_or_else(|| inconsistent("missing control heterozygous-site count"))?
            )
            .map_err(RsomicsError::Io)?;
        }
    }
    Ok(())
}

fn write_estimate(
    writer: &mut dyn Write,
    chromosome: &ChromosomeCall,
    estimate: crate::call::AberrantEstimate,
) -> Result<()> {
    let (start, end) = chromosome_bounds(chromosome)?;
    writeln!(
        writer,
        "CF\t{}\t{}\t{}\t{:.2}\t{:.6}",
        chromosome.reference_name, start, end, estimate.fraction, estimate.baf_deviation
    )
    .map_err(RsomicsError::Io)
}

fn chromosome_bounds(chromosome: &ChromosomeCall) -> Result<(u32, u32)> {
    let start = chromosome
        .sites
        .first()
        .ok_or_else(|| inconsistent("chromosome has no sites"))?
        .position;
    let end = chromosome
        .sites
        .last()
        .ok_or_else(|| inconsistent("chromosome has no sites"))?
        .position;
    Ok((start, end))
}

fn write_file(path: PathBuf, body: impl FnOnce(&mut dyn Write) -> Result<()>) -> Result<()> {
    let file = File::create(&path)
        .rs_with_context(|| format!("creating staged report {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    body(&mut writer)?;
    writer
        .flush()
        .rs_with_context(|| format!("flushing staged report {}", path.display()))?;
    writer
        .get_ref()
        .sync_all()
        .rs_with_context(|| format!("syncing staged report {}", path.display()))
}

fn validate_call_result(result: &CallResult) -> Result<()> {
    validate_sample_name(&result.sample)?;
    if let Some(control) = result.control_sample.as_deref() {
        validate_sample_name(control)?;
        if control == result.sample {
            return Err(inconsistent("query and control sample names are equal"));
        }
    }
    if result.chromosomes.is_empty() {
        return Err(inconsistent("result has no chromosomes"));
    }
    let optimized = result.chromosomes[0].query_estimate.is_some();
    let mut references = HashSet::new();
    for chromosome in &result.chromosomes {
        if chromosome.reference_name.is_empty()
            || !references.insert(chromosome.reference_name.as_str())
        {
            return Err(inconsistent("reference names are empty or duplicated"));
        }
        if chromosome.sites.is_empty() || chromosome.regions.is_empty() {
            return Err(inconsistent(format!(
                "{} has no sites or regions",
                chromosome.reference_name
            )));
        }
        if chromosome.query_estimate.is_some() != optimized {
            return Err(inconsistent(
                "query aberrant-fraction estimates are incomplete",
            ));
        }
        if result.control_sample.is_some() {
            if chromosome.control_estimate.is_some() != optimized {
                return Err(inconsistent(
                    "control aberrant-fraction estimates are incomplete",
                ));
            }
        } else if chromosome.control_estimate.is_some() {
            return Err(inconsistent(
                "control estimate is present without a control sample",
            ));
        }
        for estimate in [chromosome.query_estimate, chromosome.control_estimate]
            .into_iter()
            .flatten()
        {
            if !estimate.fraction.is_finite()
                || !(0.0..=1.0).contains(&estimate.fraction)
                || !estimate.baf_deviation.is_finite()
                || estimate.baf_deviation <= 0.0
            {
                return Err(inconsistent(format!(
                    "{} has an invalid aberrant-fraction estimate",
                    chromosome.reference_name
                )));
            }
        }
    }
    Ok(())
}

fn validate_polysomy_result(result: &PolysomyResult) -> Result<()> {
    if result.sample.is_empty() {
        return Err(inconsistent_polysomy("sample name is empty"));
    }
    if result.chromosomes.is_empty() {
        return Err(inconsistent_polysomy("result has no chromosomes"));
    }
    let mut references = HashSet::new();
    for chromosome in &result.chromosomes {
        let distribution = &chromosome.distribution;
        if distribution.reference_name.is_empty()
            || !references.insert(distribution.reference_name.as_str())
        {
            return Err(inconsistent_polysomy(
                "reference names are empty or duplicated",
            ));
        }
        if distribution.observations == 0 || distribution.bins.len() != BINS {
            return Err(inconsistent_polysomy(format!(
                "{} has an invalid observation or bin count",
                distribution.reference_name
            )));
        }
        if !(distribution.rr_boundary < distribution.heterozygous_center
            && distribution.heterozygous_center < distribution.aa_boundary
            && distribution.aa_boundary <= distribution.fitted_end
            && distribution.fitted_end < BINS)
        {
            return Err(inconsistent_polysomy(format!(
                "{} has invalid distribution boundaries",
                distribution.reference_name
            )));
        }
        for (index, bin) in distribution.bins.iter().enumerate() {
            let expected = index as f64 / (BINS - 1) as f64;
            if !bin.baf.is_finite()
                || (bin.baf - expected).abs() > 1e-12
                || !bin.normalized_count.is_finite()
                || bin.normalized_count < 0.0
            {
                return Err(inconsistent_polysomy(format!(
                    "{} has an invalid distribution bin",
                    distribution.reference_name
                )));
            }
        }
        if !chromosome.copy_number.is_finite() {
            return Err(inconsistent_polysomy(format!(
                "{} copy number is not finite",
                distribution.reference_name
            )));
        }
        validate_deviation(chromosome.absolute_deviation, &distribution.reference_name)?;
        let mut selected = None;
        for (index, candidate) in chromosome.candidates.iter().enumerate() {
            if !(2..=4).contains(&candidate.model_copy_number)
                || !candidate.estimated_copy_number.is_finite()
            {
                return Err(inconsistent_polysomy(format!(
                    "{} has an invalid candidate model",
                    distribution.reference_name
                )));
            }
            validate_deviation(candidate.absolute_deviation, &distribution.reference_name)?;
            for curve in &candidate.curves {
                if curve.start_bin > curve.end_bin
                    || curve.end_bin >= BINS
                    || !curve.absolute_deviation.is_finite()
                    || curve.absolute_deviation < 0.0
                    || curve.function.is_empty()
                    || curve.fitted.len() != curve.end_bin - curve.start_bin + 1
                {
                    return Err(inconsistent_polysomy(format!(
                        "{} has an invalid fit curve",
                        distribution.reference_name
                    )));
                }
                for (index, bin) in (curve.start_bin..=curve.end_bin).zip(&curve.fitted) {
                    let expected = index as f64 / (BINS - 1) as f64;
                    if !bin.baf.is_finite()
                        || (bin.baf - expected).abs() > 1e-12
                        || !bin.normalized_count.is_finite()
                        || bin.normalized_count < 0.0
                    {
                        return Err(inconsistent_polysomy(format!(
                            "{} has an invalid fitted curve point",
                            distribution.reference_name
                        )));
                    }
                }
            }
            if candidate.selected {
                if candidate.rejection.is_some() || selected.replace(index).is_some() {
                    return Err(inconsistent_polysomy(format!(
                        "{} has inconsistent selected candidates",
                        distribution.reference_name
                    )));
                }
            } else if candidate.rejection.is_none() {
                return Err(inconsistent_polysomy(format!(
                    "{} has an unselected candidate without a rejection",
                    distribution.reference_name
                )));
            }
        }
        match selected {
            Some(index)
                if chromosome.copy_number == chromosome.candidates[index].estimated_copy_number
                    && chromosome.absolute_deviation
                        == chromosome.candidates[index].absolute_deviation => {}
            Some(_) => {
                return Err(inconsistent_polysomy(format!(
                    "{} selected candidate differs from the final call",
                    distribution.reference_name
                )));
            }
            None if chromosome.candidates.is_empty()
                && distribution
                    .preliminary_copy_number
                    .is_some_and(|value| f64::from(value) == chromosome.copy_number)
                && chromosome.absolute_deviation.is_none() => {}
            None if !chromosome.candidates.is_empty() && chromosome.copy_number == -1.0 => {}
            None => {
                return Err(inconsistent_polysomy(format!(
                    "{} has no decision supporting the final call",
                    distribution.reference_name
                )));
            }
        }
    }
    Ok(())
}

fn validate_deviation(value: Option<f64>, reference_name: &str) -> Result<()> {
    if value.is_some_and(|value| !value.is_finite() || value < 0.0) {
        return Err(inconsistent_polysomy(format!(
            "{reference_name} has an invalid absolute deviation"
        )));
    }
    Ok(())
}

fn validate_sample_name(name: &str) -> Result<()> {
    if name.is_empty() || name == "." || name == ".." || name.contains(['/', '\\']) {
        return Err(RsomicsError::InvalidInput(format!(
            "sample name {name:?} cannot be used in report filenames"
        )));
    }
    Ok(())
}

fn inconsistent(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(format!("inconsistent call result: {}", message.into()))
}

fn inconsistent_polysomy(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(format!("inconsistent polysomy result: {}", message.into()))
}

fn scientific(value: f64) -> String {
    let raw = format!("{value:.6e}");
    let (mantissa, exponent) = raw.split_once('e').unwrap();
    let exponent = exponent.parse::<i32>().unwrap();
    format!("{mantissa}e{exponent:+03}")
}
