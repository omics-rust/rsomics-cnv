use std::path::Path;
use std::process::Command;

use rsomics_cnv::call::{CallOptions, analyze, analyze_with_allele_frequencies};
use rsomics_cnv::reports::write_call_reports;
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

fn write_optimization_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"BAF\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"LRR\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    );
    for index in 0..600 {
        let position = index * 1_000 + 1;
        let baf = [0.0, 0.4, 0.6, 1.0][index % 4];
        text.push_str(&format!(
            "chr1\t{position}\t.\tA\tG\t.\t.\t.\tBAF:LRR\t{baf:.7}:0.3000000\n"
        ));
    }
    std::fs::write(path, text).unwrap();
}

fn write_paired_optimization_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"BAF\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"LRR\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQUERY\tCONTROL\n",
    );
    for index in 0..600 {
        let position = index * 1_000 + 1;
        let query = [0.0, 0.4, 0.6, 1.0][index % 4];
        let control = [0.0, 0.357_142_857, 0.642_857_143, 1.0][index % 4];
        text.push_str(&format!(
            "chr1\t{position}\t.\tA\tG\t.\t.\t.\tBAF:LRR\t{query:.9}:0.3000000\t{control:.9}:0.3000000\n"
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
fn optimization_recovers_diluted_cn3_parameters() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("optimization.vcf");
    write_optimization_fixture(&input);
    let result = analyze(
        &input,
        SampleSelection::default(),
        CallOptions {
            optimize_aberrant_fraction: Some(0.2),
            ..CallOptions::default()
        },
    )
    .unwrap();
    let estimate = result.chromosomes[0].query_estimate.unwrap();
    assert!((estimate.fraction - 0.5).abs() < 1e-6);
    assert!((estimate.baf_deviation - 0.028_284_271).abs() < 1e-6);
    assert!(
        result.chromosomes[0]
            .sites
            .iter()
            .all(|site| site.copy_number == 3)
    );
}

#[test]
fn optimization_uses_cn3_variance_without_homozygous_alternate_sites() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("optimization.vcf");
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"BAF\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"LRR\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    );
    for index in 0..300 {
        let position = index * 1_000 + 1;
        let baf = if index % 2 == 0 { 0.4 } else { 0.6 };
        text.push_str(&format!(
            "chr1\t{position}\t.\tA\tG\t.\t.\t.\tBAF:LRR\t{baf:.7}:0.3000000\n"
        ));
    }
    std::fs::write(&input, text).unwrap();
    let result = analyze(
        &input,
        SampleSelection::default(),
        CallOptions {
            optimize_aberrant_fraction: Some(0.2),
            ..CallOptions::default()
        },
    )
    .unwrap();
    let estimate = result.chromosomes[0].query_estimate.unwrap();
    assert!((estimate.fraction - 0.5).abs() < 1e-6);
    assert!((estimate.baf_deviation - 0.028_284_271).abs() < 1e-6);
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
    let error = analyze_with_allele_frequencies(
        Path::new("missing.vcf"),
        SampleSelection::default(),
        CallOptions {
            lrr_smoothing_window: 0,
            ..CallOptions::default()
        },
        Path::new("missing.tsv"),
    )
    .unwrap_err();
    assert!(error.to_string().contains("smoothing window"), "{error}");
    let error = analyze(
        Path::new("missing.vcf"),
        SampleSelection::default(),
        CallOptions {
            optimize_aberrant_fraction: Some(f64::NAN),
            ..CallOptions::default()
        },
    )
    .unwrap_err();
    assert!(
        error.to_string().contains("optimization minimum"),
        "{error}"
    );
}

#[test]
fn allele_frequency_file_restricts_sites_and_checks_order() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("cnv.vcf");
    let frequencies = directory.path().join("frequencies.tsv");
    write_fixture(&input);
    std::fs::write(
        &frequencies,
        "#CHROM\tPOS\tREF,ALT\tAF\n\
chr1\t1000\tA,G\t0.10\n\
chr1\t6000\tA,T\t0.25\n\
chr1\t11000\tA,G\t.\n",
    )
    .unwrap();
    let result = analyze_with_allele_frequencies(
        &input,
        SampleSelection::default(),
        CallOptions::default(),
        &frequencies,
    )
    .unwrap();
    assert_eq!(result.chromosomes[0].sites.len(), 3);

    std::fs::write(
        &frequencies,
        "chr1\t6000\tA,G\t0.10\nchr1\t1000\tA,G\t0.20\n",
    )
    .unwrap();
    let error = analyze_with_allele_frequencies(
        &input,
        SampleSelection::default(),
        CallOptions::default(),
        &frequencies,
    )
    .unwrap_err();
    assert!(error.to_string().contains("coordinate sorted"), "{error}");
}

#[test]
fn result_separates_observed_and_modeled_lrr() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("smoothing.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"BAF\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"LRR\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t10\t.\tA\tG\t.\t.\t.\tBAF:LRR\t0.5:0.0\n\
chr1\t20\t.\tA\tG\t.\t.\t.\tBAF:LRR\t0.5:0.9\n\
chr1\t30\t.\tA\tG\t.\t.\t.\tBAF:LRR\t0.5:0.0\n",
    )
    .unwrap();
    let result = analyze(
        &input,
        SampleSelection::default(),
        CallOptions {
            lrr_smoothing_window: 3,
            ..CallOptions::default()
        },
    )
    .unwrap();
    let site = &result.chromosomes[0].sites[1];
    assert!((site.measurement.lrr.unwrap() - 0.9).abs() < 1e-7);
    assert!((site.modeled_lrr.unwrap() - 0.3).abs() < 1e-7);
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
    let ours_output = directory.path().join("rsomics");
    write_call_reports(&ours_output, &ours).unwrap();
    assert_eq!(
        std::fs::read(output.join("dat.SAMPLE.tab")).unwrap(),
        std::fs::read(ours_output.join("dat.SAMPLE.tab")).unwrap()
    );
    assert_eq!(
        std::fs::read(output.join("cn.SAMPLE.tab")).unwrap(),
        std::fs::read(ours_output.join("cn.SAMPLE.tab")).unwrap()
    );
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
    let ours_output = directory.path().join("rsomics");
    write_call_reports(&ours_output, &ours).unwrap();
    let chromosome = &ours.chromosomes[0];

    for (name, control) in [("QUERY", false), ("CONTROL", true)] {
        assert_eq!(
            std::fs::read(output.join(format!("dat.{name}.tab"))).unwrap(),
            std::fs::read(ours_output.join(format!("dat.{name}.tab"))).unwrap()
        );
        assert_eq!(
            std::fs::read(output.join(format!("cn.{name}.tab"))).unwrap(),
            std::fs::read(ours_output.join(format!("cn.{name}.tab"))).unwrap()
        );
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

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn allele_frequency_restriction_matches_bcftools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("cnv.vcf");
    let frequencies = directory.path().join("frequencies.tsv");
    let output = directory.path().join("bcftools");
    write_fixture(&input);
    let mut table = String::new();
    for index in (0..180).step_by(3) {
        let position = index * 5_000 + 1_000;
        let (alleles, frequency) = match index % 9 {
            0 => ("A,G", "0.25"),
            3 => ("A,T", "0.40"),
            _ => ("A,G", "."),
        };
        table.push_str(&format!("chr1\t{position}\t{alleles}\t{frequency}\n"));
    }
    std::fs::write(&frequencies, table).unwrap();
    let compression = Command::new("bgzip")
        .args(["--force"])
        .arg(&frequencies)
        .output()
        .unwrap();
    assert!(
        compression.status.success(),
        "{}",
        String::from_utf8_lossy(&compression.stderr)
    );
    let frequencies = directory.path().join("frequencies.tsv.gz");
    let indexing = Command::new("tabix")
        .args(["--sequence", "1", "--begin", "2", "--end", "2"])
        .arg(&frequencies)
        .output()
        .unwrap();
    assert!(
        indexing.status.success(),
        "{}",
        String::from_utf8_lossy(&indexing.stderr)
    );
    let upstream = Command::new("bcftools")
        .args(["cnv", "-s", "SAMPLE", "-f"])
        .arg(&frequencies)
        .arg("-o")
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        upstream.status.success(),
        "{}",
        String::from_utf8_lossy(&upstream.stderr)
    );
    let ours = analyze_with_allele_frequencies(
        &input,
        SampleSelection {
            query: Some("SAMPLE".to_owned()),
            control: None,
        },
        CallOptions::default(),
        &frequencies,
    )
    .unwrap();
    let ours_output = directory.path().join("rsomics");
    write_call_reports(&ours_output, &ours).unwrap();
    assert_eq!(
        std::fs::read(output.join("dat.SAMPLE.tab")).unwrap(),
        std::fs::read(ours_output.join("dat.SAMPLE.tab")).unwrap()
    );
    assert_eq!(
        std::fs::read(output.join("cn.SAMPLE.tab")).unwrap(),
        std::fs::read(ours_output.join("cn.SAMPLE.tab")).unwrap()
    );
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn aberrant_fraction_optimization_matches_bcftools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("optimization.vcf");
    let upstream_output = directory.path().join("bcftools");
    write_optimization_fixture(&input);
    let upstream = Command::new("bcftools")
        .args(["cnv", "-s", "SAMPLE", "-O", "0.2", "-o"])
        .arg(&upstream_output)
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
        SampleSelection::default(),
        CallOptions {
            optimize_aberrant_fraction: Some(0.2),
            ..CallOptions::default()
        },
    )
    .unwrap();
    let ours_output = directory.path().join("rsomics");
    write_call_reports(&ours_output, &ours).unwrap();
    assert_eq!(
        std::fs::read(upstream_output.join("dat.SAMPLE.tab")).unwrap(),
        std::fs::read(ours_output.join("dat.SAMPLE.tab")).unwrap()
    );
    assert_eq!(
        std::fs::read(upstream_output.join("cn.SAMPLE.tab")).unwrap(),
        std::fs::read(ours_output.join("cn.SAMPLE.tab")).unwrap()
    );
    let cf = |path: &Path| {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .find(|line| line.starts_with("CF\t"))
            .unwrap()
            .to_owned()
    };
    assert_eq!(
        cf(&upstream_output.join("summary.SAMPLE.tab")),
        cf(&ours_output.join("summary.SAMPLE.tab"))
    );
    let estimate = ours.chromosomes[0].query_estimate.unwrap();
    assert!((estimate.fraction - 0.5).abs() < 1e-6);
    assert!((estimate.baf_deviation - 0.028284).abs() < 1e-6);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn paired_optimization_matches_bcftools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("optimization.vcf");
    let upstream_output = directory.path().join("bcftools");
    write_paired_optimization_fixture(&input);
    let upstream = Command::new("bcftools")
        .args(["cnv", "-s", "QUERY", "-c", "CONTROL", "-O", "0.2", "-o"])
        .arg(&upstream_output)
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
        CallOptions {
            optimize_aberrant_fraction: Some(0.2),
            ..CallOptions::default()
        },
    )
    .unwrap();
    let ours_output = directory.path().join("rsomics");
    write_call_reports(&ours_output, &ours).unwrap();
    for name in ["QUERY", "CONTROL"] {
        assert_eq!(
            std::fs::read(upstream_output.join(format!("cn.{name}.tab"))).unwrap(),
            std::fs::read(ours_output.join(format!("cn.{name}.tab"))).unwrap()
        );
        let cf = |root: &Path| {
            std::fs::read_to_string(root.join(format!("summary.{name}.tab")))
                .unwrap()
                .lines()
                .find(|line| line.starts_with("CF\t"))
                .unwrap()
                .to_owned()
        };
        assert_eq!(cf(&upstream_output), cf(&ours_output));
    }
    let upstream_joint = std::fs::read_to_string(upstream_output.join("summary.tab")).unwrap();
    let ours_joint = std::fs::read_to_string(ours_output.join("summary.tab")).unwrap();
    assert_eq!(
        upstream_joint.lines().find(|line| line.starts_with("CF\t")),
        ours_joint.lines().find(|line| line.starts_with("CF\t"))
    );
}
