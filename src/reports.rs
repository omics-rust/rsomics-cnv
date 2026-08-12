use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError};
use serde::Serialize;

use crate::call::{CallResult, ChromosomeCall};

const DAT_HEADER: &str = "# [1]Chromosome\t[2]Position\t[3]BAF\t[4]LRR\n";
const CN_HEADER: &str =
    "# [1]Chromosome\t[2]Position\t[3]CN\t[4]P(CN0)\t[5]P(CN1)\t[6]P(CN2)\t[7]P(CN3)\n";
const SUMMARY_HEADER: &str = "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy Number state\t[6]Quality\t[7]nSites\t[8]nHETs\n";

#[derive(Serialize)]
struct ResultDocument<'a> {
    schema: &'static str,
    producer: Producer,
    result: &'a CallResult,
}

#[derive(Serialize)]
struct Producer {
    name: &'static str,
    version: &'static str,
}

pub fn write_call_reports(output: &Path, result: &CallResult) -> Result<()> {
    validate_result(result)?;
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

    write_sample_reports(stage.path(), result, false)?;
    if result.control_sample.is_some() {
        write_sample_reports(stage.path(), result, true)?;
        write_file(stage.path().join("summary.tab"), |writer| {
            write_joint_summary(writer, result)
        })?;
    }
    write_file(stage.path().join("result.json"), |writer| {
        let document = ResultDocument {
            schema: "rsomics-cnv/call-result/v1",
            producer: Producer {
                name: "rsomics-cnv",
                version: env!("CARGO_PKG_VERSION"),
            },
            result,
        };
        serde_json::to_writer_pretty(&mut *writer, &document)
            .map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
        writer.write_all(b"\n").map_err(RsomicsError::Io)
    })?;
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
    for chromosome in chromosomes {
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
    for chromosome in &result.chromosomes {
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

fn validate_result(result: &CallResult) -> Result<()> {
    validate_sample_name(&result.sample)?;
    if let Some(control) = result.control_sample.as_deref() {
        validate_sample_name(control)?;
        if control == result.sample {
            return Err(inconsistent("query and control sample names are equal"));
        }
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
