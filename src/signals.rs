use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;
use noodles_bcf as bcf;
use noodles_vcf::{
    self as vcf,
    header::record::value::map::format,
    variant::{RecordBuf, record_buf::samples::sample::Value as SampleValue},
};
use rsomics_common::{Context, Result, RsomicsError};

type Source = Box<dyn Read>;
type Buffered = BufReader<Source>;
type Compressed = MultiGzDecoder<Buffered>;

enum Input {
    Bcf(bcf::io::Reader<Compressed>),
    BcfRaw(bcf::io::Reader<Buffered>),
    Vcf(Buffered),
    VcfGz(BufReader<Compressed>),
}

struct VariantReader {
    input: Input,
    source: String,
}

#[derive(Default)]
struct RecordScratch {
    bcf: bcf::Record,
    text: Vec<u8>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RequiredSignals {
    Baf,
    BafAndLrr,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SampleSelection {
    pub query: Option<String>,
    pub control: Option<String>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Measurement {
    pub baf: Option<f64>,
    pub lrr: Option<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SignalSite {
    pub reference_name: String,
    pub position: u32,
    pub query: Measurement,
    pub control: Option<Measurement>,
}

pub struct SignalReader {
    reader: VariantReader,
    header: vcf::Header,
    scratch: RecordScratch,
    required: RequiredSignals,
    query: usize,
    control: Option<usize>,
    query_name: String,
    control_name: Option<String>,
    previous: Option<(usize, u32)>,
    records: u64,
}

impl fmt::Debug for SignalReader {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SignalReader")
            .field("source", &self.reader.source)
            .field("required", &self.required)
            .field("query_name", &self.query_name)
            .field("control_name", &self.control_name)
            .field("records", &self.records)
            .finish()
    }
}

impl SignalReader {
    pub fn open(path: &Path, samples: SampleSelection, required: RequiredSignals) -> Result<Self> {
        let mut reader = VariantReader::open(path)?;
        let header = reader.read_header()?;
        validate_format(&header, "BAF")?;
        if required == RequiredSignals::BafAndLrr {
            validate_format(&header, "LRR")?;
        }

        let (query, query_name) = select_query(&header, samples.query.as_deref())?;
        let (control, control_name) = match samples.control.as_deref() {
            Some(name) => {
                let index = select_named(&header, name, "control")?;
                if index == query {
                    return Err(invalid("query and control samples must be different"));
                }
                (Some(index), Some(name.to_owned()))
            }
            None => (None, None),
        };

        Ok(Self {
            reader,
            header,
            scratch: RecordScratch::default(),
            required,
            query,
            control,
            query_name,
            control_name,
            previous: None,
            records: 0,
        })
    }

    pub fn query_sample(&self) -> &str {
        &self.query_name
    }

    pub fn control_sample(&self) -> Option<&str> {
        self.control_name.as_deref()
    }

    pub fn next_site(&mut self) -> Result<Option<SignalSite>> {
        loop {
            let number = self.records + 1;
            let Some(record) = self
                .reader
                .read_record(&self.header, &mut self.scratch, number)?
            else {
                return Ok(None);
            };
            self.records = number;

            let reference_name = record.reference_sequence_name().to_owned();
            let reference = self
                .header
                .contigs()
                .get_index_of(&reference_name)
                .ok_or_else(|| {
                    record_error(
                        &self.reader.source,
                        number,
                        format!("reference {reference_name:?} is absent from the header"),
                    )
                })?;
            let position = record
                .variant_start()
                .map(usize::from)
                .and_then(|value| u32::try_from(value).ok())
                .ok_or_else(|| {
                    record_error(
                        &self.reader.source,
                        number,
                        "record position is absent or exceeds uint32",
                    )
                })?;
            self.check_order(reference, position, number)?;

            let query = measurement(
                &record,
                self.query,
                self.required,
                &self.reader.source,
                number,
                &self.query_name,
            )?;
            let control = self
                .control
                .map(|index| {
                    measurement(
                        &record,
                        index,
                        self.required,
                        &self.reader.source,
                        number,
                        self.control_name.as_deref().unwrap(),
                    )
                })
                .transpose()?;
            if query.baf.is_none() && control.is_none_or(|value| value.baf.is_none()) {
                continue;
            }

            return Ok(Some(SignalSite {
                reference_name,
                position,
                query,
                control,
            }));
        }
    }

    fn check_order(&mut self, reference: usize, position: u32, number: u64) -> Result<()> {
        if self
            .previous
            .is_some_and(|(previous_reference, previous_position)| {
                reference < previous_reference
                    || (reference == previous_reference && position < previous_position)
            })
        {
            return Err(record_error(
                &self.reader.source,
                number,
                "records must be coordinate sorted in header contig order",
            ));
        }
        self.previous = Some((reference, position));
        Ok(())
    }
}

fn validate_format(header: &vcf::Header, name: &str) -> Result<()> {
    let definition = header
        .formats()
        .get(name)
        .ok_or_else(|| invalid(format!("FORMAT/{name} is absent from the header")))?;
    if definition.number() != format::Number::Count(1) || definition.ty() != format::Type::Float {
        return Err(invalid(format!(
            "FORMAT/{name} must have Number=1,Type=Float"
        )));
    }
    Ok(())
}

fn select_query(header: &vcf::Header, name: Option<&str>) -> Result<(usize, String)> {
    if let Some(name) = name {
        return select_named(header, name, "query").map(|index| (index, name.to_owned()));
    }
    match header.sample_names().len() {
        0 => Err(invalid("VCF or BCF has no samples")),
        1 => Ok((0, header.sample_names().get_index(0).unwrap().to_owned())),
        count => Err(invalid(format!(
            "input has {count} samples; choose a query sample"
        ))),
    }
}

fn select_named(header: &vcf::Header, name: &str, role: &str) -> Result<usize> {
    header
        .sample_names()
        .get_index_of(name)
        .ok_or_else(|| invalid(format!("unknown {role} sample {name:?}")))
}

fn measurement(
    record: &RecordBuf,
    sample: usize,
    required: RequiredSignals,
    source: &str,
    number: u64,
    sample_name: &str,
) -> Result<Measurement> {
    let baf = sample_float(record, sample, "BAF", source, number, sample_name)?;
    if let Some(value) = baf
        && !(0.0..=1.0).contains(&value)
    {
        return Err(record_error(
            source,
            number,
            format!("sample {sample_name:?} FORMAT/BAF must be between 0 and 1"),
        ));
    }

    match required {
        RequiredSignals::Baf => Ok(Measurement { baf, lrr: None }),
        RequiredSignals::BafAndLrr => {
            let lrr = sample_float(record, sample, "LRR", source, number, sample_name)?;
            if baf.is_some() && lrr.is_some() {
                Ok(Measurement { baf, lrr })
            } else {
                Ok(Measurement::default())
            }
        }
    }
}

fn sample_float(
    record: &RecordBuf,
    sample: usize,
    name: &str,
    source: &str,
    number: u64,
    sample_name: &str,
) -> Result<Option<f64>> {
    let Some(values) = record.samples().select(name) else {
        return Ok(None);
    };
    let Some(value) = values.get(sample).flatten() else {
        return Ok(None);
    };
    let SampleValue::Float(value) = value else {
        return Err(record_error(
            source,
            number,
            format!("sample {sample_name:?} FORMAT/{name} is not a float"),
        ));
    };
    let value = f64::from(*value);
    if !value.is_finite() {
        return Err(record_error(
            source,
            number,
            format!("sample {sample_name:?} must have a finite {name} value"),
        ));
    }
    Ok(Some(value))
}

impl VariantReader {
    fn open(path: &Path) -> Result<Self> {
        let source = if path == Path::new("-") {
            "standard input".to_owned()
        } else {
            path.display().to_string()
        };
        let input: Source = if path == Path::new("-") {
            Box::new(io::stdin())
        } else {
            Box::new(
                File::open(path)
                    .rs_with_context(|| format!("opening variant input {}", path.display()))?,
            )
        };
        let mut input = BufReader::new(input);
        let compressed = is_gzip(&mut input)
            .map_err(|error| input_error(&source, "detecting compression", error))?;
        let bcf = is_bcf(&mut input, compressed)
            .map_err(|error| input_error(&source, "detecting variant format", error))?;
        let input = match (bcf, compressed) {
            (true, true) => Input::Bcf(bcf::io::Reader::from(MultiGzDecoder::new(input))),
            (true, false) => Input::BcfRaw(bcf::io::Reader::from(input)),
            (false, true) => Input::VcfGz(BufReader::new(MultiGzDecoder::new(input))),
            (false, false) => Input::Vcf(input),
        };
        Ok(Self { input, source })
    }

    fn read_header(&mut self) -> Result<vcf::Header> {
        let raw = match &mut self.input {
            Input::Bcf(reader) => read_bcf_header(reader),
            Input::BcfRaw(reader) => read_bcf_header(reader),
            Input::Vcf(reader) => read_vcf_header(reader),
            Input::VcfGz(reader) => read_vcf_header(reader),
        }
        .map_err(|error| input_error(&self.source, "reading VCF header", error))?;
        let text = std::str::from_utf8(&raw).map_err(|error| {
            input_error(
                &self.source,
                "decoding VCF header",
                io::Error::new(io::ErrorKind::InvalidData, error),
            )
        })?;
        let mut header: vcf::Header = text.parse().map_err(|error| {
            input_error(
                &self.source,
                "parsing VCF header",
                io::Error::new(io::ErrorKind::InvalidData, error),
            )
        })?;
        if !matches!(self.input, Input::Vcf(_) | Input::VcfGz(_)) {
            *header.string_maps_mut() = text.parse().map_err(|error| {
                input_error(
                    &self.source,
                    "reading BCF string maps",
                    io::Error::new(io::ErrorKind::InvalidData, error),
                )
            })?;
        }
        Ok(header)
    }

    fn read_record(
        &mut self,
        header: &vcf::Header,
        scratch: &mut RecordScratch,
        number: u64,
    ) -> Result<Option<RecordBuf>> {
        let index =
            usize::try_from(number).map_err(|_| invalid("variant record count exceeds usize"))?;
        match &mut self.input {
            Input::Vcf(reader) => {
                read_vcf_record(reader, header, &mut scratch.text, number, &self.source)
            }
            Input::VcfGz(reader) => {
                read_vcf_record(reader, header, &mut scratch.text, number, &self.source)
            }
            Input::Bcf(reader) => read_bcf_record(
                reader,
                header,
                &mut scratch.bcf,
                index,
                number,
                &self.source,
            ),
            Input::BcfRaw(reader) => read_bcf_record(
                reader,
                header,
                &mut scratch.bcf,
                index,
                number,
                &self.source,
            ),
        }
    }
}

fn read_vcf_record(
    reader: &mut impl BufRead,
    header: &vcf::Header,
    buffer: &mut Vec<u8>,
    number: u64,
    source: &str,
) -> Result<Option<RecordBuf>> {
    buffer.clear();
    if reader
        .read_until(b'\n', buffer)
        .map_err(|error| input_error(source, &format!("reading variant record {number}"), error))?
        == 0
    {
        return Ok(None);
    }
    trim_line_ending(buffer);
    let record = vcf::Record::try_from(buffer.as_slice())
        .map_err(|error| record_error(source, number, format!("parsing VCF: {error}")))?;
    RecordBuf::try_from_variant_record(header, &record)
        .map(Some)
        .map_err(|error| record_error(source, number, format!("decoding VCF: {error}")))
}

fn read_bcf_record<R>(
    reader: &mut bcf::io::Reader<R>,
    header: &vcf::Header,
    record: &mut bcf::Record,
    index: usize,
    number: u64,
    source: &str,
) -> Result<Option<RecordBuf>>
where
    R: Read,
{
    if reader
        .read_record(record)
        .map_err(|error| input_error(source, &format!("reading variant record {index}"), error))?
        == 0
    {
        return Ok(None);
    }
    RecordBuf::try_from_variant_record(header, record)
        .map(Some)
        .map_err(|error| record_error(source, number, format!("decoding BCF: {error}")))
}

fn is_gzip(reader: &mut impl BufRead) -> io::Result<bool> {
    Ok(reader
        .fill_buf()?
        .get(..2)
        .is_some_and(|magic| magic == [0x1f, 0x8b]))
}

fn is_bcf(reader: &mut impl BufRead, compressed: bool) -> io::Result<bool> {
    let source = reader.fill_buf()?;
    if compressed {
        let mut decoder = MultiGzDecoder::new(source);
        let mut magic = [0; 3];
        decoder.read_exact(&mut magic)?;
        Ok(magic == *b"BCF")
    } else {
        Ok(source.get(..3).is_some_and(|magic| magic == b"BCF"))
    }
}

fn read_vcf_header(reader: &mut impl BufRead) -> io::Result<Vec<u8>> {
    let mut raw = Vec::new();
    loop {
        let source = reader.fill_buf()?;
        if source.first() != Some(&b'#') {
            break;
        }
        reader.read_until(b'\n', &mut raw)?;
    }
    Ok(raw)
}

fn read_bcf_header<R>(reader: &mut bcf::io::Reader<R>) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let mut reader = reader.header_reader();
    if reader.read_magic_number()? != *b"BCF" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BCF magic number",
        ));
    }
    let version = reader.read_format_version()?;
    if version != (2, 2) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported BCF version {}.{}", version.0, version.1),
        ));
    }
    let mut raw_reader = reader.raw_vcf_header_reader()?;
    let mut raw = Vec::new();
    raw_reader.read_to_end(&mut raw)?;
    raw_reader.discard_to_end()?;
    Ok(raw)
}

fn trim_line_ending(line: &mut Vec<u8>) {
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

fn input_error(source: &str, action: &str, error: io::Error) -> RsomicsError {
    invalid(format!("{source}: {action}: {error}"))
}

fn record_error(source: &str, number: u64, message: impl fmt::Display) -> RsomicsError {
    invalid(format!("{source}: record {number}: {message}"))
}
