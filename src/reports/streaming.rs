use super::*;

struct ReportFile {
    path: PathBuf,
    writer: BufWriter<File>,
}

struct SampleReportFiles {
    data: ReportFile,
    copy_number: ReportFile,
    summary: ReportFile,
}

pub(crate) struct CallReportWriter {
    stage: tempfile::TempDir,
    output: PathBuf,
    sample: String,
    control_sample: Option<String>,
    query: SampleReportFiles,
    control: Option<SampleReportFiles>,
    joint_summary: Option<ReportFile>,
    chromosomes: Vec<CallChromosomeReport>,
    references: HashSet<String>,
    optimized: Option<bool>,
    options: CallReportOptions,
}

impl ReportFile {
    fn create(path: PathBuf) -> Result<Self> {
        let file = File::create(&path)
            .rs_with_context(|| format!("creating staged report {}", path.display()))?;
        Ok(Self {
            path,
            writer: BufWriter::new(file),
        })
    }

    fn writer(&mut self) -> &mut dyn Write {
        &mut self.writer
    }

    fn finish(mut self) -> Result<()> {
        self.writer
            .flush()
            .rs_with_context(|| format!("flushing staged report {}", self.path.display()))?;
        self.writer
            .get_ref()
            .sync_all()
            .rs_with_context(|| format!("syncing staged report {}", self.path.display()))
    }
}

impl SampleReportFiles {
    fn create(directory: &Path, sample: &str) -> Result<Self> {
        let mut data = ReportFile::create(directory.join(format!("dat.{sample}.tab")))?;
        data.writer()
            .write_all(DAT_HEADER.as_bytes())
            .map_err(RsomicsError::Io)?;
        let mut copy_number = ReportFile::create(directory.join(format!("cn.{sample}.tab")))?;
        copy_number
            .writer()
            .write_all(CN_HEADER.as_bytes())
            .map_err(RsomicsError::Io)?;
        Ok(Self {
            data,
            copy_number,
            summary: ReportFile::create(directory.join(format!("summary.{sample}.tab")))?,
        })
    }

    fn finish(self) -> Result<()> {
        self.data.finish()?;
        self.copy_number.finish()?;
        self.summary.finish()
    }
}

impl CallReportWriter {
    pub(crate) fn new(
        output: &Path,
        sample: &str,
        control_sample: Option<&str>,
        options: CallReportOptions,
    ) -> Result<Self> {
        validate_sample_name(sample)?;
        if let Some(control) = control_sample {
            validate_sample_name(control)?;
            if control == sample {
                return Err(inconsistent("query and control sample names are equal"));
            }
        }
        validate_plot_threshold(options.plot_threshold)?;
        let stage = create_report_directory(output)?;
        let query = SampleReportFiles::create(stage.path(), sample)?;
        let control = control_sample
            .map(|control| SampleReportFiles::create(stage.path(), control))
            .transpose()?;
        let joint_summary = control_sample
            .map(|_| ReportFile::create(stage.path().join("summary.tab")))
            .transpose()?;
        Ok(Self {
            stage,
            output: output.to_owned(),
            sample: sample.to_owned(),
            control_sample: control_sample.map(str::to_owned),
            query,
            control,
            joint_summary,
            chromosomes: Vec::new(),
            references: HashSet::new(),
            optimized: None,
            options,
        })
    }

    pub(crate) fn write_chromosome(&mut self, chromosome: &ChromosomeCall) -> Result<()> {
        let optimized = chromosome.query_estimate.is_some();
        validate_call_chromosome(
            chromosome,
            self.control_sample.is_some(),
            self.optimized,
            &mut self.references,
        )?;
        if self.optimized.is_none() {
            self.write_summary_headers(optimized)?;
            self.optimized = Some(optimized);
        }
        write_data_chromosome(self.query.data.writer(), chromosome, false)?;
        write_copy_number_chromosome(self.query.copy_number.writer(), chromosome, false)?;
        write_summary_chromosome(self.query.summary.writer(), chromosome, false)?;
        if let Some(control) = &mut self.control {
            write_data_chromosome(control.data.writer(), chromosome, true)?;
            write_copy_number_chromosome(control.copy_number.writer(), chromosome, true)?;
            write_summary_chromosome(control.summary.writer(), chromosome, true)?;
        }
        if let Some(joint) = &mut self.joint_summary {
            write_joint_summary_chromosome(joint.writer(), chromosome)?;
        }
        if let Some(threshold) = self.options.plot_threshold
            && let Some((name, svg)) = call_plot(
                chromosome,
                &self.sample,
                self.control_sample.as_deref(),
                threshold,
            )?
        {
            write_plot(self.stage.path(), name, svg)?;
        }
        self.chromosomes.push(CallChromosomeReport {
            reference_name: chromosome.reference_name.clone(),
            sites: chromosome.sites.len(),
            regions: chromosome.regions.clone(),
            query_estimate: chromosome.query_estimate,
            control_estimate: chromosome.control_estimate,
        });
        Ok(())
    }

    fn write_summary_headers(&mut self, optimized: bool) -> Result<()> {
        write_sample_summary_header(self.query.summary.writer(), optimized)?;
        if let Some(control) = &mut self.control {
            write_sample_summary_header(control.summary.writer(), optimized)?;
        }
        if let Some(joint) = &mut self.joint_summary {
            write_joint_summary_header(
                joint.writer(),
                &self.sample,
                self.control_sample
                    .as_deref()
                    .ok_or_else(|| inconsistent("missing control sample"))?,
                optimized,
            )?;
        }
        Ok(())
    }

    pub(crate) fn finish(self) -> Result<()> {
        if self.chromosomes.is_empty() {
            return Err(inconsistent("result has no chromosomes"));
        }
        let Self {
            stage,
            output,
            sample,
            control_sample,
            query,
            control,
            joint_summary,
            chromosomes,
            options,
            ..
        } = self;
        query.finish()?;
        if let Some(control) = control {
            control.finish()?;
        }
        if let Some(joint) = joint_summary {
            joint.finish()?;
        }
        let sample_artifacts = |sample: &str| SampleArtifacts {
            data: format!("dat.{sample}.tab"),
            copy_number: format!("cn.{sample}.tab"),
            summary: format!("summary.{sample}.tab"),
        };
        let report = CallReport {
            artifacts: CallArtifacts {
                query: sample_artifacts(&sample),
                control: control_sample.as_deref().map(sample_artifacts),
                joint_summary: control_sample.as_ref().map(|_| "summary.tab".to_owned()),
            },
            sample,
            control_sample,
            chromosomes,
            plot_threshold: options.plot_threshold,
        };
        write_json(
            stage.path().join("result.json"),
            "rsomics-cnv/call-result/v2",
            &report,
        )?;
        commit_report_directory(stage, &output)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::SiteCall;
    use crate::signals::Measurement;

    fn chromosome(reference_name: &str) -> ChromosomeCall {
        ChromosomeCall {
            reference_name: reference_name.to_owned(),
            sites: vec![SiteCall {
                position: 11,
                copy_number: 2,
                control_copy_number: Some(1),
                state_probability: 0.94,
                posterior: [0.01, 0.02, 0.94, 0.03],
                control_posterior: Some([0.01, 0.91, 0.06, 0.02]),
                measurement: Measurement {
                    baf: Some(0.52),
                    lrr: Some(-0.12),
                },
                modeled_lrr: Some(-0.08),
                control_measurement: Some(Measurement {
                    baf: Some(0.48),
                    lrr: Some(-0.31),
                }),
                control_modeled_lrr: Some(-0.27),
            }],
            regions: vec![RegionCall {
                start: 11,
                end: 11,
                copy_number: 2,
                control_copy_number: Some(1),
                quality: 10.4,
                sites: 1,
                heterozygous_sites: 1,
                control_heterozygous_sites: Some(1),
            }],
            query_estimate: None,
            control_estimate: None,
        }
    }

    #[test]
    fn commits_complete_bundle() {
        let directory = tempfile::tempdir().unwrap();
        let output = directory.path().join("reports");
        let mut writer = CallReportWriter::new(
            &output,
            "QUERY",
            Some("CONTROL"),
            CallReportOptions::default(),
        )
        .unwrap();
        writer.write_chromosome(&chromosome("chr1")).unwrap();
        writer.finish().unwrap();

        assert!(output.join("dat.QUERY.tab").is_file());
        assert!(output.join("cn.CONTROL.tab").is_file());
        assert!(output.join("summary.tab").is_file());
        let json: serde_json::Value =
            serde_json::from_slice(&fs::read(output.join("result.json")).unwrap()).unwrap();
        assert_eq!(json["schema"], "rsomics-cnv/call-result/v2");
        assert_eq!(json["result"]["chromosomes"][0]["reference_name"], "chr1");
        assert_eq!(json["result"]["chromosomes"][0]["sites"], 1);
    }

    #[test]
    fn drop_leaves_no_destination() {
        let directory = tempfile::tempdir().unwrap();
        let output = directory.path().join("reports");
        let writer =
            CallReportWriter::new(&output, "QUERY", None, CallReportOptions::default()).unwrap();
        drop(writer);
        assert!(!output.exists());
    }
}
