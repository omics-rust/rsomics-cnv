use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Lines};
use std::path::{Path, PathBuf};

use flate2::bufread::MultiGzDecoder;
use rsomics_common::{Context, Result, RsomicsError};

use crate::emission::AlleleFrequencies;

const DEFAULT_NON_REFERENCE_FREQUENCY: f64 = 0.1;

struct Entry {
    reference: usize,
    position: u32,
    alleles: Vec<String>,
    non_reference_frequency: Option<f64>,
}

pub(crate) struct AlleleFrequencyReader {
    lines: Lines<Box<dyn BufRead>>,
    path: PathBuf,
    references: HashMap<String, usize>,
    current: Option<Entry>,
    previous: Option<(usize, u32)>,
    line: usize,
}

impl AlleleFrequencyReader {
    pub(crate) fn open(path: &Path, references: &[String]) -> Result<Self> {
        let file = File::open(path)
            .rs_with_context(|| format!("opening allele-frequency file {}", path.display()))?;
        let mut source = BufReader::new(file);
        let compressed = source
            .fill_buf()
            .rs_with_context(|| format!("reading allele-frequency file {}", path.display()))?
            .starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn BufRead> = if compressed {
            Box::new(BufReader::new(MultiGzDecoder::new(source)))
        } else {
            Box::new(source)
        };
        let mut frequencies = Self {
            lines: reader.lines(),
            path: path.to_path_buf(),
            references: references
                .iter()
                .enumerate()
                .map(|(index, name)| (name.clone(), index))
                .collect(),
            current: None,
            previous: None,
            line: 0,
        };
        frequencies.advance()?;
        if frequencies.current.is_none() {
            return Err(RsomicsError::InvalidInput(format!(
                "allele-frequency file {} has no records",
                path.display()
            )));
        }
        Ok(frequencies)
    }

    pub(crate) fn frequencies(
        &mut self,
        reference: usize,
        position: u32,
        alleles: &[String],
    ) -> Result<Option<AlleleFrequencies>> {
        while self
            .current
            .as_ref()
            .is_some_and(|entry| (entry.reference, entry.position) < (reference, position))
        {
            self.advance()?;
        }
        let Some(entry) = self
            .current
            .as_ref()
            .filter(|entry| (entry.reference, entry.position) == (reference, position))
        else {
            return Ok(None);
        };
        let frequency = if entry.alleles == alleles {
            entry
                .non_reference_frequency
                .unwrap_or(DEFAULT_NON_REFERENCE_FREQUENCY)
        } else {
            DEFAULT_NON_REFERENCE_FREQUENCY
        };
        AlleleFrequencies::from_non_reference_frequency(frequency).map(Some)
    }

    pub(crate) fn finish(mut self) -> Result<()> {
        while self.current.is_some() {
            self.advance()?;
        }
        Ok(())
    }

    fn advance(&mut self) -> Result<()> {
        self.current = None;
        for line in self.lines.by_ref() {
            self.line += 1;
            let line = line.rs_with_context(|| {
                format!(
                    "reading allele-frequency file {} line {}",
                    self.path.display(),
                    self.line
                )
            })?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() != 4 {
                return Err(self.record_error(
                    "expected CHROM, POS, REF,ALT, and AF in four tab-separated columns",
                ));
            }
            let reference = self.references.get(fields[0]).copied().ok_or_else(|| {
                self.record_error(format!(
                    "reference {:?} is absent from the variant header",
                    fields[0]
                ))
            })?;
            let position = fields[1]
                .parse::<u32>()
                .map_err(|_| self.record_error(format!("invalid position {:?}", fields[1])))?;
            if position == 0 {
                return Err(self.record_error("position must be positive"));
            }
            if self
                .previous
                .is_some_and(|previous| (reference, position) <= previous)
            {
                return Err(self.record_error(
                    "records must be coordinate sorted in variant-header order with unique positions",
                ));
            }
            self.previous = Some((reference, position));
            let alleles = fields[2].split(',').map(str::to_owned).collect::<Vec<_>>();
            if alleles.len() < 2 || alleles.iter().any(String::is_empty) {
                return Err(self.record_error(
                    "allele column must contain comma-separated REF and ALT values",
                ));
            }
            let non_reference_frequency = if fields[3] == "." {
                None
            } else {
                let value = fields[3].parse::<f64>().map_err(|_| {
                    self.record_error(format!("invalid allele frequency {:?}", fields[3]))
                })?;
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    return Err(
                        self.record_error("allele frequency must be finite and between 0 and 1")
                    );
                }
                Some(value)
            };
            self.current = Some(Entry {
                reference,
                position,
                alleles,
                non_reference_frequency,
            });
            return Ok(());
        }
        Ok(())
    }

    fn record_error(&self, message: impl Into<String>) -> RsomicsError {
        RsomicsError::InvalidInput(format!(
            "{}:{}: {}",
            self.path.display(),
            self.line,
            message.into()
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::{Compression, write::GzEncoder};

    use super::*;

    #[test]
    fn compressed_input_streams_in_header_order() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("frequencies.tsv.gz");
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder
            .write_all(b"chr1\t10\tA,G\t0.25\nchr1\t20\tC,T\t.\n")
            .unwrap();
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();
        let mut reader = AlleleFrequencyReader::open(&path, &["chr1".to_owned()]).unwrap();
        let exact = reader
            .frequencies(0, 10, &["A".to_owned(), "G".to_owned()])
            .unwrap()
            .unwrap();
        assert_eq!(
            exact,
            AlleleFrequencies::new(0.5625, 0.375, 0.0625).unwrap()
        );
        let mismatch = reader
            .frequencies(0, 20, &["C".to_owned(), "G".to_owned()])
            .unwrap()
            .unwrap();
        assert_eq!(
            mismatch,
            AlleleFrequencies::from_non_reference_frequency(0.1).unwrap()
        );
        assert!(reader.frequencies(0, 30, &[]).unwrap().is_none());
        reader.finish().unwrap();
    }
}
