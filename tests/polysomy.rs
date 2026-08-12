use std::path::Path;
use std::process::Command;

use rsomics_cnv::polysomy::{PolysomyOptions, analyze_distributions};

fn write_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
##contig=<ID=chr2,length=1000000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B-allele frequency\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    );
    for index in 0..500 {
        let baf = [0.0, 0.5, 1.0][index % 3];
        text.push_str(&format!(
            "chr1\t{}\t.\tA\tG\t.\t.\t.\tBAF\t{baf:.7}\n",
            index + 1
        ));
    }
    for index in 0..600 {
        let baf = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0][index % 4];
        text.push_str(&format!(
            "chr2\t{}\t.\tA\tG\t.\t.\t.\tBAF\t{baf:.7}\n",
            index + 1
        ));
    }
    std::fs::write(path, text).unwrap();
}

#[test]
fn distributions_are_typed_and_chromosome_partitioned() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("polysomy.vcf");
    write_fixture(&input);
    let result = analyze_distributions(&input, None, PolysomyOptions::default()).unwrap();
    assert_eq!(result.sample, "SAMPLE");
    assert_eq!(result.chromosomes.len(), 2);
    assert_eq!(result.chromosomes[0].reference_name, "chr1");
    assert_eq!(result.chromosomes[0].observations, 500);
    assert_eq!(result.chromosomes[1].reference_name, "chr2");
    assert_eq!(result.chromosomes[1].observations, 600);
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
        .args(["polysomy", "-s", "SAMPLE", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        upstream.status.success(),
        "{}",
        String::from_utf8_lossy(&upstream.stderr)
    );
    let result = analyze_distributions(
        &input,
        Some("SAMPLE".to_owned()),
        PolysomyOptions::default(),
    )
    .unwrap();
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
}
