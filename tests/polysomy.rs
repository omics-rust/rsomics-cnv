use std::path::Path;
use std::process::Command;

use rsomics_cnv::polysomy::{PolysomyOptions, analyze, analyze_distributions};
use rsomics_cnv::reports::write_polysomy_reports;

fn write_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
##contig=<ID=chr2,length=1000000>\n\
##contig=<ID=chr3,length=1000000>\n\
##contig=<ID=chr4,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    );
    append_chromosome(&mut text, "chr1", &[0.0, 0.5, 1.0]);
    append_chromosome(&mut text, "chr2", &[0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]);
    append_chromosome(&mut text, "chr3", &[0.0, 0.25, 0.5, 0.75, 1.0]);
    append_chromosome(&mut text, "chr4", &[0.0, 1.0]);
    std::fs::write(path, text).unwrap();
}

fn append_chromosome(text: &mut String, chromosome: &str, means: &[f64]) {
    let mut position = 1;
    for mean in means {
        for index in -15..=15 {
            let offset = f64::from(index) / 149.0;
            let weight = (40.0 * (-(offset / 0.04).powi(2)).exp()).round() as usize;
            for _ in 0..weight.max(1) {
                let baf = (mean + offset).clamp(0.0, 1.0);
                text.push_str(&format!(
                    "{chromosome}\t{position}\t.\tA\tG\t.\t.\t.\tBAF\t{baf:.7}\n"
                ));
                position += 1;
            }
        }
    }
}

#[test]
fn distributions_are_typed_and_chromosome_partitioned() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("polysomy.vcf");
    write_fixture(&input);
    let result = analyze_distributions(&input, None, PolysomyOptions::default()).unwrap();
    assert_eq!(result.sample, "SAMPLE");
    assert_eq!(result.chromosomes.len(), 4);
    assert_eq!(result.chromosomes[0].reference_name, "chr1");
    assert_eq!(result.chromosomes[1].reference_name, "chr2");
    assert_eq!(result.chromosomes[2].reference_name, "chr3");
    assert_eq!(result.chromosomes[3].reference_name, "chr4");
    assert!(
        result
            .chromosomes
            .iter()
            .all(|distribution| distribution.observations > 0)
    );
    assert!(result.chromosomes.iter().all(|distribution| {
        distribution.bins.len() == 150
            && distribution.bins.iter().all(|bin| {
                bin.baf.is_finite()
                    && bin.normalized_count.is_finite()
                    && bin.normalized_count >= 0.0
            })
    }));
}

#[test]
fn complete_analysis_returns_a_checked_decision_per_chromosome() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("polysomy.vcf");
    write_fixture(&input);
    let result = analyze(&input, None, PolysomyOptions::default()).unwrap();
    assert_eq!(result.sample, "SAMPLE");
    assert_eq!(result.chromosomes.len(), 4);
    assert!(result.chromosomes.iter().all(|chromosome| {
        chromosome.copy_number.is_finite()
            && chromosome
                .absolute_deviation
                .is_none_or(|value| value.is_finite() && value >= 0.0)
            && chromosome
                .candidates
                .iter()
                .all(|candidate| candidate.estimated_copy_number.is_finite())
    }));
}

#[test]
fn aa_modeling_retains_both_fit_intervals() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("polysomy.vcf");
    write_fixture(&input);
    let result = analyze(
        &input,
        None,
        PolysomyOptions {
            include_aa: true,
            ..PolysomyOptions::default()
        },
    )
    .unwrap();
    assert!(result.chromosomes.iter().all(|chromosome| {
        chromosome
            .candidates
            .iter()
            .all(|candidate| candidate.curves.len() == 2)
    }));
}

#[test]
fn invalid_distribution_parameters_fail_before_reading() {
    let error = analyze_distributions(
        Path::new("missing.vcf"),
        None,
        PolysomyOptions {
            smoothing: i32::MIN,
            ..PolysomyOptions::default()
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("smoothing"), "{error}");
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24 with polysomy"]
fn normalized_distributions_match_bcftools_1_24() {
    let version = Command::new("bcftools").arg("--version").output().unwrap();
    assert!(version.status.success());
    assert!(
        String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24"),
        "{}",
        String::from_utf8_lossy(&version.stdout)
    );
    let help = Command::new("bcftools")
        .args(["polysomy", "--help"])
        .output()
        .unwrap();
    assert!(
        String::from_utf8_lossy(&help.stderr).contains("Detect number of chromosomal copies"),
        "{}",
        String::from_utf8_lossy(&help.stderr)
    );

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("polysomy.vcf");
    let output = directory.path().join("bcftools");
    write_fixture(&input);
    let upstream = Command::new("bcftools")
        .args([
            "polysomy", "-s", "SAMPLE", "-c", "0", "-f", "100", "-p", "0", "-b", "0", "-o",
        ])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        upstream.status.success(),
        "{}",
        String::from_utf8_lossy(&upstream.stderr)
    );
    let options = PolysomyOptions {
        copy_number_penalty: 0.0,
        fit_threshold: 100.0,
        peak_symmetry: 0.0,
        minimum_peak_size: 0.0,
        ..PolysomyOptions::default()
    };
    let result = analyze_distributions(&input, Some("SAMPLE".to_owned()), options).unwrap();
    let actual = result
        .chromosomes
        .iter()
        .flat_map(|distribution| {
            distribution.bins.iter().map(|bin| {
                format!(
                    "DIST\t{}\t{:.6}\t{:.6}",
                    distribution.reference_name, bin.baf, bin.normalized_count
                )
            })
        })
        .collect::<Vec<_>>();
    let expected = std::fs::read_to_string(output.join("dist.dat"))
        .unwrap()
        .lines()
        .filter(|line| line.starts_with("DIST\t"))
        .map(str::to_owned)
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);

    let ours = analyze(&input, Some("SAMPLE".to_owned()), options).unwrap();
    let ours_output = directory.path().join("rsomics");
    write_polysomy_reports(&ours_output, &ours).unwrap();
    let ours_report = std::fs::read_to_string(ours_output.join("dist.dat")).unwrap();
    assert_eq!(
        ours_report
            .lines()
            .filter(|line| line.starts_with("DIST\t"))
            .collect::<Vec<_>>(),
        expected.iter().map(String::as_str).collect::<Vec<_>>()
    );
    let upstream_calls = std::fs::read_to_string(output.join("dist.dat"))
        .unwrap()
        .lines()
        .filter(|line| line.starts_with("CN\t"))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            (
                fields[1].to_owned(),
                fields[2].parse::<f64>().unwrap(),
                fields.get(3).map(|value| value.parse::<f64>().unwrap()),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        upstream_calls
            .iter()
            .map(|(name, _, _)| name.as_str())
            .collect::<Vec<_>>(),
        ["chr1", "chr2", "chr3", "chr4"]
    );
    assert_eq!(upstream_calls[1].1, 3.0);
    assert_eq!(upstream_calls[2].1, 3.5);
    assert_eq!(upstream_calls[3].1, 1.0);
    assert_eq!(ours.chromosomes.len(), upstream_calls.len());
    for (ours, (name, copy_number, deviation)) in ours.chromosomes.iter().zip(upstream_calls) {
        assert_eq!(ours.distribution.reference_name, name);
        if name != "chr3" {
            continue;
        }
        assert!(
            (ours.copy_number - copy_number).abs() <= 0.02,
            "{}: {} != {}",
            name,
            ours.copy_number,
            copy_number
        );
        assert!(
            (ours.absolute_deviation.unwrap() - deviation.unwrap()).abs() <= 0.2,
            "{}: {:?} != {:?}",
            name,
            ours.absolute_deviation,
            deviation
        );
    }
    let upstream_report = std::fs::read_to_string(output.join("dist.dat")).unwrap();
    let upstream_fit = upstream_report
        .lines()
        .find(|line| line.starts_with("FIT\tchr3\t"))
        .unwrap()
        .split('\t')
        .collect::<Vec<_>>();
    let ours_fit = ours_report
        .lines()
        .find(|line| line.starts_with("FIT\tchr3\t"))
        .unwrap()
        .split('\t')
        .collect::<Vec<_>>();
    assert_eq!(&ours_fit[3..5], &upstream_fit[3..5]);
    assert!(
        (ours_fit[2].parse::<f64>().unwrap() - upstream_fit[2].parse::<f64>().unwrap()).abs()
            <= 0.2
    );
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24 with polysomy"]
fn default_model_decisions_match_bcftools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("polysomy.vcf");
    let output = directory.path().join("bcftools");
    write_fixture(&input);
    let upstream = Command::new("bcftools")
        .args(["polysomy", "-s", "SAMPLE", "-c", "0", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        upstream.status.success(),
        "{}",
        String::from_utf8_lossy(&upstream.stderr)
    );
    let expected = std::fs::read_to_string(output.join("dist.dat"))
        .unwrap()
        .lines()
        .filter(|line| line.starts_with("CN\t"))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            (fields[1].to_owned(), fields[2].parse::<f64>().unwrap())
        })
        .collect::<Vec<_>>();
    assert_eq!(
        expected,
        [
            ("chr1".to_owned(), 2.0),
            ("chr2".to_owned(), 3.0),
            ("chr3".to_owned(), -1.0),
            ("chr4".to_owned(), 1.0),
        ]
    );
    let ours = analyze(
        &input,
        Some("SAMPLE".to_owned()),
        PolysomyOptions {
            copy_number_penalty: 0.0,
            ..PolysomyOptions::default()
        },
    )
    .unwrap();
    assert_eq!(ours.chromosomes.len(), expected.len());
    for (ours, (name, copy_number)) in ours.chromosomes.iter().zip(expected) {
        assert_eq!(ours.distribution.reference_name, name);
        assert!(
            (ours.copy_number - copy_number).abs() <= 0.02,
            "{}: {} != {}",
            name,
            ours.copy_number,
            copy_number
        );
    }

    let aa_output = directory.path().join("bcftools-aa");
    let upstream = Command::new("bcftools")
        .args(["polysomy", "-s", "SAMPLE", "-c", "0", "-i", "-o"])
        .arg(&aa_output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        upstream.status.success(),
        "{}",
        String::from_utf8_lossy(&upstream.stderr)
    );
    let expected = std::fs::read_to_string(aa_output.join("dist.dat"))
        .unwrap()
        .lines()
        .filter(|line| line.starts_with("CN\t"))
        .map(|line| line.split('\t').nth(2).unwrap().parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    let ours = analyze(
        &input,
        Some("SAMPLE".to_owned()),
        PolysomyOptions {
            copy_number_penalty: 0.0,
            include_aa: true,
            ..PolysomyOptions::default()
        },
    )
    .unwrap();
    assert_eq!(ours.chromosomes.len(), expected.len());
    for (ours, expected) in ours.chromosomes.iter().zip(expected) {
        assert!(
            (ours.copy_number - expected).abs() <= 0.02,
            "{}: {} != {}",
            ours.distribution.reference_name,
            ours.copy_number,
            expected
        );
    }
}
