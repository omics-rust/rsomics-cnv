use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_cnv::signals::{Measurement, RequiredSignals, SampleSelection, SignalReader};

const VCF: &str = "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=1000>\n\
##contig=<ID=chr2,length=1000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R ratio\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\tnormal\n\
chr1\t10\t.\tA\tG\t.\tPASS\t.\tBAF:LRR\t0.25:-0.1\t0.50:0.0\n\
chr1\t20\t.\tC\tT\t.\tPASS\t.\tBAF:LRR\t.:.\t0.75:0.2\n\
chr2\t5\t.\tG\tA\t.\tPASS\t.\tBAF:LRR\t1.0:0.3\t.:.\n";

fn write_representations(directory: &Path) -> [PathBuf; 3] {
    let vcf_path = directory.join("signals.vcf");
    let bgzf_path = directory.join("signals.vcf.gz");
    let bcf_path = directory.join("signals.bcf");
    std::fs::write(&vcf_path, VCF).unwrap();

    let mut bgzf_writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
    bgzf_writer.write_all(VCF.as_bytes()).unwrap();
    bgzf_writer.try_finish().unwrap();

    let mut reader = vcf::io::Reader::new(BufReader::new(File::open(&vcf_path).unwrap()));
    let header = reader.read_header().unwrap();
    let mut bcf_writer = bcf::io::Writer::new(File::create(&bcf_path).unwrap());
    bcf_writer.write_header(&header).unwrap();
    for record in reader.records() {
        bcf_writer
            .write_variant_record(&header, &record.unwrap())
            .unwrap();
    }
    bcf_writer.try_finish().unwrap();

    [vcf_path, bgzf_path, bcf_path]
}

fn selection() -> SampleSelection {
    SampleSelection {
        query: Some("tumor".to_owned()),
        control: Some("normal".to_owned()),
    }
}

#[test]
fn all_variant_encodings_share_the_signal_contract() {
    let directory = tempfile::tempdir().unwrap();
    let paths = write_representations(directory.path());
    let mut expected = None;

    for path in paths {
        let mut reader =
            SignalReader::open(&path, selection(), RequiredSignals::BafAndLrr).unwrap();
        assert_eq!(reader.query_sample(), "tumor");
        assert_eq!(reader.control_sample(), Some("normal"));

        let mut sites = Vec::new();
        while let Some(site) = reader.next_site().unwrap() {
            sites.push(site);
        }
        if let Some(expected) = &expected {
            assert_eq!(&sites, expected, "{}", path.display());
        } else {
            expected = Some(sites);
        }
    }

    let sites = expected.unwrap();
    assert_eq!(sites.len(), 3);
    assert_eq!(sites[0].reference_name, "chr1");
    assert_eq!(sites[0].position, 10);
    assert_eq!(
        sites[0].query,
        Measurement {
            baf: Some(0.25),
            lrr: Some(f64::from(-0.1_f32)),
        }
    );
    assert_eq!(sites[1].query, Measurement::default());
    assert_eq!(sites[2].control, Some(Measurement::default()));
}

#[test]
fn sample_and_schema_boundaries_fail_loudly() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("signals.vcf");
    std::fs::write(&input, VCF).unwrap();

    let error =
        SignalReader::open(&input, SampleSelection::default(), RequiredSignals::Baf).unwrap_err();
    assert!(
        error.to_string().contains("choose a query sample"),
        "{error}"
    );

    let error = SignalReader::open(
        &input,
        SampleSelection {
            query: Some("missing".to_owned()),
            control: None,
        },
        RequiredSignals::Baf,
    )
    .unwrap_err();
    assert!(
        error.to_string().contains("unknown query sample"),
        "{error}"
    );

    let invalid = directory.path().join("invalid-schema.vcf");
    std::fs::write(
        &invalid,
        VCF.replace("ID=BAF,Number=1,Type=Float", "ID=BAF,Number=2,Type=Integer"),
    )
    .unwrap();
    let error = SignalReader::open(&invalid, selection(), RequiredSignals::Baf).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("FORMAT/BAF must have Number=1,Type=Float"),
        "{error}"
    );
}

#[test]
fn invalid_values_and_record_order_are_rejected_with_context() {
    let directory = tempfile::tempdir().unwrap();
    let out_of_range = directory.path().join("out-of-range.vcf");
    std::fs::write(&out_of_range, VCF.replace("0.25:-0.1", "1.25:-0.1")).unwrap();
    let mut reader =
        SignalReader::open(&out_of_range, selection(), RequiredSignals::BafAndLrr).unwrap();
    let error = reader.next_site().unwrap_err();
    assert!(error.to_string().contains("record 1"), "{error}");
    assert!(error.to_string().contains("BAF"), "{error}");

    let non_finite = directory.path().join("non-finite.vcf");
    std::fs::write(&non_finite, VCF.replace("0.25:-0.1", "0.25:nan")).unwrap();
    let mut reader =
        SignalReader::open(&non_finite, selection(), RequiredSignals::BafAndLrr).unwrap();
    let error = reader.next_site().unwrap_err();
    assert!(error.to_string().contains("finite LRR"), "{error}");

    let unsorted = directory.path().join("unsorted.vcf");
    std::fs::write(&unsorted, VCF.replace("chr1\t20", "chr1\t5")).unwrap();
    let mut reader =
        SignalReader::open(&unsorted, selection(), RequiredSignals::BafAndLrr).unwrap();
    reader.next_site().unwrap().unwrap();
    let error = reader.next_site().unwrap_err();
    assert!(error.to_string().contains("coordinate sorted"), "{error}");
}

#[test]
fn baf_only_mode_does_not_require_lrr() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("baf.vcf");
    let baf_only = VCF
        .replace(
            "##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R ratio\">\n",
            "",
        )
        .replace("BAF:LRR", "BAF")
        .replace("0.25:-0.1", "0.25")
        .replace("0.50:0.0", "0.50")
        .replace(".:.", ".")
        .replace("0.75:0.2", "0.75")
        .replace("1.0:0.3", "1.0");
    std::fs::write(&input, baf_only).unwrap();

    let mut reader = SignalReader::open(&input, selection(), RequiredSignals::Baf).unwrap();
    let site = reader.next_site().unwrap().unwrap();
    assert_eq!(site.query.lrr, None);
}
