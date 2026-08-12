use std::path::Path;
use std::process::Command;

use rsomics_cnv::call::{CallOptions, analyze};
use rsomics_cnv::signals::SampleSelection;

fn write_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R ratio\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    );
    for index in 0..180 {
        let position = index * 5_000 + 1_000;
        let (baf, lrr) = match index / 60 {
            0 => ([0.0, 0.5, 1.0][index % 3], 0.0),
            1 => ([0.0, 1.0][index % 2], -0.45),
            _ => ([0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0][index % 4], 0.3),
        };
        text.push_str(&format!(
            "chr1\t{position}\t.\tA\tG\t.\tPASS\t.\tBAF:LRR\t{baf:.7}:{lrr:.7}\n"
        ));
    }
    std::fs::write(path, text).unwrap();
}

fn write_paired_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R ratio\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQUERY\tCONTROL\n",
    );
    for index in 0..120 {
        let position = index * 5_000 + 1_000;
        let (query_baf, query_lrr) = if index < 60 {
            ([0.0, 0.5, 1.0][index % 3], 0.0)
        } else {
            ([0.0, 1.0][index % 2], -0.45)
        };
        let control_baf = [0.0, 0.5, 1.0][index % 3];
        text.push_str(&format!(
            "chr1\t{position}\t.\tA\tG\t.\tPASS\t.\tBAF:LRR\t{query_baf:.7}:{query_lrr:.7}\t{control_baf:.7}:0.0000000\n"
        ));
    }
    std::fs::write(path, text).unwrap();
}

#[test]
fn synthetic_chromosome_recovers_three_copy_number_states() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("cnv.vcf");
    write_fixture(&input);

    let result = analyze(
        &input,
        SampleSelection {
            query: Some("SAMPLE".to_owned()),
            control: None,
        },
        CallOptions::default(),
    )
    .unwrap();
    assert_eq!(result.sample, "SAMPLE");
    assert_eq!(result.chromosomes.len(), 1);
    assert_eq!(result.chromosomes[0].sites.len(), 180);
    assert_eq!(
        result.chromosomes[0]
            .regions
            .iter()
            .map(|region| region.copy_number)
            .collect::<Vec<_>>(),
        [2, 1, 3]
    );
    assert!(
        result.chromosomes[0]
            .regions
            .iter()
            .all(|region| region.quality.is_finite())
    );
}

#[test]
fn invalid_call_parameters_fail_before_reading() {
    let error = analyze(
        Path::new("missing.vcf"),
        SampleSelection::default(),
        CallOptions {
            lrr_smoothing_window: 0,
            ..CallOptions::default()
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("smoothing window"), "{error}");
}

#[test]
fn paired_analysis_retains_query_and_control_marginals() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("paired.vcf");
    write_paired_fixture(&input);
    let result = analyze(
        &input,
        SampleSelection {
            query: Some("QUERY".to_owned()),
            control: Some("CONTROL".to_owned()),
        },
        CallOptions::default(),
    )
    .unwrap();
    assert_eq!(result.control_sample.as_deref(), Some("CONTROL"));
    let chromosome = &result.chromosomes[0];
    assert_eq!(
        chromosome
            .regions
            .iter()
            .map(|region| (region.copy_number, region.control_copy_number))
            .collect::<Vec<_>>(),
        [(2, Some(2)), (1, Some(2))]
    );
    assert!(chromosome.sites.iter().all(|site| {
        site.control_posterior.is_some()
            && site.control_copy_number == Some(2)
            && site.state_probability.is_finite()
    }));
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn site_and_region_results_match_bcftools_1_24() {
    let version = Command::new("bcftools").arg("--version").output().unwrap();
    assert!(version.status.success());
    assert!(
        String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24"),
        "{}",
        String::from_utf8_lossy(&version.stdout)
    );

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("cnv.vcf");
    let output = directory.path().join("bcftools");
    write_fixture(&input);

    let upstream = Command::new("bcftools")
        .args(["cnv", "-s", "SAMPLE", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        upstream.status.success(),
        "{}",
        String::from_utf8_lossy(&upstream.stderr)
    );

    let ours = analyze(
        &input,
        SampleSelection {
            query: Some("SAMPLE".to_owned()),
            control: None,
        },
        CallOptions::default(),
    )
    .unwrap();
    let chromosome = &ours.chromosomes[0];
    let cn = std::fs::read_to_string(output.join("cn.SAMPLE.tab")).unwrap();
    let upstream_sites = cn
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(chromosome.sites.len(), upstream_sites.len());
    for (ours, upstream) in chromosome.sites.iter().zip(upstream_sites) {
        let fields = upstream.split('\t').collect::<Vec<_>>();
        assert_eq!(ours.position.to_string(), fields[1]);
        assert_eq!(ours.copy_number.to_string(), fields[2]);
        for (actual, expected) in ours.posterior.iter().zip(&fields[3..]) {
            let expected = expected.parse::<f64>().unwrap();
            assert!((actual - expected).abs() < 1e-6, "{actual} != {expected}");
        }
    }

    let summary = std::fs::read_to_string(output.join("summary.SAMPLE.tab")).unwrap();
    let upstream_regions = summary
        .lines()
        .filter(|line| line.starts_with("RG\t"))
        .collect::<Vec<_>>();
    assert_eq!(chromosome.regions.len(), upstream_regions.len());
    for (ours, upstream) in chromosome.regions.iter().zip(upstream_regions) {
        let fields = upstream.split('\t').collect::<Vec<_>>();
        assert_eq!(ours.start.to_string(), fields[2]);
        assert_eq!(ours.end.to_string(), fields[3]);
        assert_eq!(ours.copy_number.to_string(), fields[4]);
        assert_eq!(
            ours.sites,
            chromosome
                .sites
                .iter()
                .filter(|site| {
                    site.position >= ours.start
                        && site.position <= ours.end
                        && site.copy_number == ours.copy_number
                })
                .count()
        );
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn paired_marginals_match_bcftools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("paired.vcf");
    let output = directory.path().join("bcftools");
    write_paired_fixture(&input);
    let upstream = Command::new("bcftools")
        .args(["cnv", "-s", "QUERY", "-c", "CONTROL", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        upstream.status.success(),
        "{}",
        String::from_utf8_lossy(&upstream.stderr)
    );
    let ours = analyze(
        &input,
        SampleSelection {
            query: Some("QUERY".to_owned()),
            control: Some("CONTROL".to_owned()),
        },
        CallOptions::default(),
    )
    .unwrap();
    let chromosome = &ours.chromosomes[0];

    for (name, control) in [("QUERY", false), ("CONTROL", true)] {
        let cn = std::fs::read_to_string(output.join(format!("cn.{name}.tab"))).unwrap();
        let upstream = cn
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect::<Vec<_>>();
        assert_eq!(chromosome.sites.len(), upstream.len());
        for (ours, upstream) in chromosome.sites.iter().zip(upstream) {
            let fields = upstream.split('\t').collect::<Vec<_>>();
            let state = if control {
                ours.control_copy_number.unwrap()
            } else {
                ours.copy_number
            };
            let posterior = if control {
                ours.control_posterior.unwrap()
            } else {
                ours.posterior
            };
            assert_eq!(state.to_string(), fields[2]);
            for (actual, expected) in posterior.iter().zip(&fields[3..]) {
                let expected = expected.parse::<f64>().unwrap();
                assert!((actual - expected).abs() < 1e-6, "{actual} != {expected}");
            }
        }
    }
}
