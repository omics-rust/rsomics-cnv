use std::path::Path;
use std::process::Command;

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-cnv"))
}

fn write_call_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R ratio\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    );
    for index in 0..90 {
        let (baf, lrr) = match index / 30 {
            0 => ([0.0, 0.5, 1.0][index % 3], 0.0),
            1 => ([0.0, 1.0][index % 2], -0.45),
            _ => ([0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0][index % 4], 0.3),
        };
        text.push_str(&format!(
            "chr1\t{}\t.\tA\tG\t.\t.\t.\tBAF:LRR\t{baf:.7}:{lrr:.7}\n",
            index * 5_000 + 1_000
        ));
    }
    std::fs::write(path, text).unwrap();
}

fn write_polysomy_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    );
    let mut position = 1;
    for mean in [0.0, 0.5, 1.0] {
        for index in -15..=15 {
            let offset = f64::from(index) / 149.0;
            let weight = (40.0 * (-(offset / 0.04).powi(2)).exp()).round() as usize;
            for _ in 0..weight.max(1) {
                let baf = (mean + offset).clamp(0.0, 1.0);
                text.push_str(&format!(
                    "chr1\t{position}\t.\tA\tG\t.\t.\t.\tBAF\t{baf:.7}\n"
                ));
                position += 1;
            }
        }
    }
    std::fs::write(path, text).unwrap();
}

#[test]
fn help_exposes_one_product_tree() {
    let output = binary().arg("--help").output().unwrap();
    assert!(output.status.success());
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(help.contains("call"), "{help}");
    assert!(help.contains("polysomy"), "{help}");
    assert!(help.contains("--json"), "{help}");

    let output = binary().args(["help", "call"]).output().unwrap();
    assert!(output.status.success());
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(help.contains("--sample <NAME>"), "{help}");
    assert!(help.contains("--control <NAME>"), "{help}");
    assert!(help.contains("--allele-frequencies <TSV>"), "{help}");
    assert!(help.contains("--regions-file <FILE>"), "{help}");
    assert!(help.contains("--targets-file <FILE>"), "{help}");
    assert!(help.contains("--optimize <FRACTION>"), "{help}");
    assert!(
        help.contains("--regions-overlap <REGIONS_OVERLAP>"),
        "{help}"
    );
    assert!(help.contains("--output <DIRECTORY>"), "{help}");

    let output = binary().args(["help", "polysomy"]).output().unwrap();
    assert!(output.status.success());
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(help.contains("--regions <REGIONS>"), "{help}");
    assert!(help.contains("--targets <TARGETS>"), "{help}");
}

#[test]
fn call_writes_reports_and_shared_json_envelope() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("call.vcf");
    let frequencies = directory.path().join("frequencies.tsv");
    let targets = directory.path().join("targets.txt");
    let reports = directory.path().join("call-reports");
    write_call_fixture(&input);
    std::fs::write(
        &frequencies,
        "chr1\t1000\tA,G\t0.1\nchr1\t6000\tA,G\t0.2\nchr1\t11000\tA,G\t0.3\n",
    )
    .unwrap();
    std::fs::write(&targets, "chr1\t6000\n").unwrap();
    let output = binary()
        .args([
            "--json",
            "call",
            "--sample",
            "SAMPLE",
            "--allele-frequencies",
        ])
        .arg(&frequencies)
        .arg("--targets-file")
        .arg(&targets)
        .arg("--output")
        .arg(&reports)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let document: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(document["status"], "ok");
    assert_eq!(document["tool"], "rsomics-cnv");
    assert_eq!(document["result"]["workflow"], "call");
    assert_eq!(document["result"]["sample"], "SAMPLE");
    assert_eq!(document["result"]["sites"], 1);
    assert!(reports.join("result.json").is_file());
    assert!(reports.join("summary.SAMPLE.tab").is_file());
}

#[test]
fn polysomy_writes_reports_and_fails_on_existing_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("polysomy.vcf");
    let reports = directory.path().join("polysomy-reports");
    write_polysomy_fixture(&input);
    let run = || {
        binary()
            .args(["polysomy", "--sample", "SAMPLE", "--output"])
            .arg(&reports)
            .arg(&input)
            .output()
            .unwrap()
    };
    let output = run();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
    assert!(reports.join("dist.dat").is_file());
    assert!(reports.join("result.json").is_file());

    std::fs::remove_file(&input).unwrap();
    let output = run();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("already exists"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
