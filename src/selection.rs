use std::path::{Path, PathBuf};

use noodles_core::{Position, Region, region::Interval};
use noodles_vcf::variant::{RecordBuf, record::info::field::key, record_buf::info::field::Value};
use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::{GenomicRegion, RegionFileMode, read_region_file};

/// How a VCF or BCF record is compared with a selected interval.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OverlapMode {
    /// Match the record POS only.
    Position,
    /// Match the POS-to-END or reference-allele span.
    Record,
    /// Match the reference span changed by a literal variant.
    Variant,
}

/// Inline and file-backed site restrictions shared by CNV workflows.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SiteSelection {
    /// Inline indexed regions.
    pub regions: Vec<String>,
    /// BED, VCF, or tabular indexed regions.
    pub regions_file: Option<PathBuf>,
    /// Indexed-region overlap rule.
    pub regions_overlap: OverlapMode,
    /// Inline streaming targets.
    pub targets: Vec<String>,
    /// BED, VCF, or tabular streaming targets.
    pub targets_file: Option<PathBuf>,
    /// Streaming-target overlap rule.
    pub targets_overlap: OverlapMode,
    /// Exclude records matching targets instead of including them.
    pub exclude_targets: bool,
}

impl Default for SiteSelection {
    fn default() -> Self {
        Self {
            regions: Vec::new(),
            regions_file: None,
            regions_overlap: OverlapMode::Record,
            targets: Vec::new(),
            targets_file: None,
            targets_overlap: OverlapMode::Position,
            exclude_targets: false,
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct CompiledSelection {
    regions: Option<RegionSet>,
    targets: Option<RegionSet>,
    exclude_targets: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct RegionSet {
    regions: Vec<Region>,
    overlap: OverlapMode,
}

impl CompiledSelection {
    pub(crate) fn new(selection: &SiteSelection) -> Result<Self> {
        let regions = compile(
            &selection.regions,
            selection.regions_file.as_deref(),
            RegionFileMode::Intervals,
            selection.regions_overlap,
            "regions",
        )?;
        let targets = compile(
            &selection.targets,
            selection.targets_file.as_deref(),
            RegionFileMode::Targets,
            selection.targets_overlap,
            "targets",
        )?;
        if selection.exclude_targets && targets.is_none() {
            return Err(invalid("target exclusion requires a target list"));
        }
        Ok(Self {
            regions,
            targets,
            exclude_targets: selection.exclude_targets,
        })
    }

    pub(crate) fn keeps(&self, record: &RecordBuf) -> bool {
        self.regions
            .as_ref()
            .is_none_or(|regions| regions.matches(record))
            && self
                .targets
                .as_ref()
                .is_none_or(|targets| targets.matches(record) != self.exclude_targets)
    }

    pub(crate) fn regions(&self) -> Option<&RegionSet> {
        self.regions.as_ref()
    }
}

impl RegionSet {
    fn new(regions: Vec<Region>, overlap: OverlapMode, kind: &str) -> Result<Self> {
        if regions.is_empty() {
            return Err(invalid(format!("{kind} list is empty")));
        }
        Ok(Self { regions, overlap })
    }

    pub(crate) fn matches(&self, record: &RecordBuf) -> bool {
        self.regions.iter().any(|region| {
            let name: &[u8] = region.name().as_ref();
            name == record.reference_sequence_name().as_bytes()
                && overlaps(record, region.interval(), self.overlap)
        })
    }

    pub(crate) fn merged(&self, header: &noodles_vcf::Header) -> Result<Vec<Region>> {
        let mut values = Vec::with_capacity(self.regions.len());
        for region in &self.regions {
            let name = std::str::from_utf8(region.name().as_ref())
                .map_err(|_| invalid(format!("region reference name is not UTF-8: {region}")))?;
            let reference = header.contigs().get_index_of(name).ok_or_else(|| {
                invalid(format!(
                    "region reference is absent from the VCF header: {name}"
                ))
            })?;
            values.push((
                reference,
                name.to_owned(),
                region.interval().start().unwrap_or(Position::MIN),
                region.interval().end().unwrap_or(Position::MAX),
            ));
        }
        values.sort_by_key(|(reference, _, start, end)| (*reference, *start, *end));

        let mut merged: Vec<(usize, String, Position, Position)> = Vec::new();
        for (reference, name, start, end) in values {
            if let Some((last_reference, _, _, last_end)) = merged.last_mut()
                && *last_reference == reference
                && start <= last_end.checked_add(1).unwrap_or(Position::MAX)
            {
                *last_end = (*last_end).max(end);
                continue;
            }
            merged.push((reference, name, start, end));
        }
        Ok(merged
            .into_iter()
            .map(|(_, name, start, end)| Region::new(name, start..=end))
            .collect())
    }

    pub(crate) fn overlap(&self) -> OverlapMode {
        self.overlap
    }
}

fn compile(
    inline: &[String],
    file: Option<&Path>,
    file_mode: RegionFileMode,
    overlap: OverlapMode,
    kind: &str,
) -> Result<Option<RegionSet>> {
    if !inline.is_empty() && file.is_some() {
        return Err(invalid(format!(
            "inline {kind} and a {kind} file cannot be used together"
        )));
    }
    let regions = if let Some(path) = file {
        read_region_file(path, file_mode)
            .map_err(|error| invalid(error.to_string()))?
            .into_iter()
            .map(convert_region)
            .collect::<Result<Vec<_>>>()?
    } else if !inline.is_empty() {
        inline
            .iter()
            .flat_map(|value| value.split(','))
            .map(parse_region)
            .collect::<Result<Vec<_>>>()?
    } else {
        return Ok(None);
    };
    RegionSet::new(regions, overlap, kind).map(Some)
}

fn convert_region(region: GenomicRegion) -> Result<Region> {
    let name = region.chrom().to_owned();
    let (start, end) = match (region.start(), region.end()) {
        (None, None) => return Ok(Region::new(name, ..)),
        (Some(start), Some(end)) => (start, end),
        _ => unreachable!("genomic region bounds are paired"),
    };
    let start = usize::try_from(start)
        .ok()
        .and_then(|value| value.checked_add(1))
        .and_then(Position::new)
        .ok_or_else(|| invalid("region start exceeds the supported coordinate range"))?;
    let end = usize::try_from(end)
        .ok()
        .and_then(Position::new)
        .ok_or_else(|| invalid("region end exceeds the supported coordinate range"))?;
    Ok(Region::new(name, start..=end))
}

fn parse_region(value: &str) -> Result<Region> {
    let interval = value.rsplit_once(':').map(|(_, interval)| interval);
    let open_end = interval
        .and_then(|interval| interval.strip_suffix('-'))
        .is_some_and(is_coordinate);
    let normalized = if open_end {
        value.strip_suffix('-').expect("open-ended region")
    } else {
        value
    };
    let region: Region = normalized
        .parse()
        .map_err(|error| invalid(format!("invalid region {value:?}: {error}")))?;
    if interval.is_some_and(is_coordinate) {
        let start = region
            .interval()
            .start()
            .ok_or_else(|| invalid(format!("invalid region {value:?}: missing position")))?;
        Ok(Region::new(region.name().to_owned(), start..=start))
    } else {
        Ok(region)
    }
}

fn is_coordinate(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || byte == b',')
}

pub(crate) fn overlaps(record: &RecordBuf, interval: Interval, mode: OverlapMode) -> bool {
    let Some(start) = record.variant_start() else {
        return false;
    };
    match mode {
        OverlapMode::Position => interval.contains(start),
        OverlapMode::Record => interval.intersects(record_interval(record, start)),
        OverlapMode::Variant => {
            variant_interval(record, start).is_some_and(|variant| interval.intersects(variant))
        }
    }
}

fn record_interval(record: &RecordBuf, start: Position) -> Interval {
    let end = record
        .info()
        .get(key::END_POSITION)
        .flatten()
        .and_then(|value| match value {
            Value::Integer(value) => usize::try_from(*value).ok(),
            _ => None,
        })
        .and_then(|value| Position::try_from(value).ok())
        .filter(|end| *end >= start)
        .unwrap_or_else(|| {
            start
                .checked_add(record.reference_bases().len().max(1) - 1)
                .unwrap_or(Position::MAX)
        });
    Interval::from(start..=end)
}

fn variant_interval(record: &RecordBuf, start: Position) -> Option<Interval> {
    let alternates = record.alternate_bases();
    if alternates.as_ref().is_empty() {
        return None;
    }
    let reference = record.reference_bases().as_bytes();
    let end = record_interval(record, start)
        .end()
        .expect("bounded record");
    let mut offset = usize::from(end) - usize::from(start) + 1;
    for alternate in alternates.as_ref() {
        let prefix = reference
            .iter()
            .zip(alternate.as_bytes())
            .take_while(|(reference, alternate)| reference == alternate)
            .count();
        offset = offset.min(prefix);
        if offset == 0 {
            break;
        }
    }
    let variant_start = start.checked_add(offset)?;
    (variant_start <= end).then(|| Interval::from(variant_start..=end))
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
