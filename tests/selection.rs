use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::process::Command;

use rsomics_cnv::call::{CallOptions, analyze_selected, analyze_with_allele_frequencies_selected};
use rsomics_cnv::polysomy::{
    PolysomyOptions, analyze_distributions_selected, analyze_selected as analyze_polysomy_selected,
};
use rsomics_cnv::reports::{write_call_reports, write_polysomy_reports};
use rsomics_cnv::selection::{OverlapMode, SiteSelection};
use rsomics_cnv::signals::SampleSelection;

fn write_fixture(path: &Path) {
    std::fs::write(
        path,
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"BAF\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"LRR\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t10\t.\tAACGT\tAATGT\t.\t.\t.\tBAF:LRR\t0.5:0.0\n\
chr1\t20\t.\tA\tG\t.\t.\t.\tBAF:LRR\t0.5:0.0\n",
    )
    .unwrap();
}

fn write_indexed_vcf(plain: &Path, output: &Path) {
    let mut writer = noodles_bgzf::io::Writer::new(File::create(output).unwrap());
    writer.write_all(&std::fs::read(plain).unwrap()).unwrap();
    writer.finish().unwrap();
    let index = noodles_vcf::fs::index(output).unwrap();
    noodles_tabix::fs::write(format!("{}.tbi", output.display()), &index).unwrap();
}

fn write_polysomy_fixture(path: &Path) {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000>\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"BAF\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
    );
    for position in 1..=300 {
        let baf = [0.0, 0.5, 1.0][position % 3];
        text.push_str(&format!("chr1\t{position}\t.\tA\tG\t.\t.\t.\tBAF\t{baf}\n"));
    }
    std::fs::write(path, text).unwrap();
}

fn assert_cn_close(upstream: &Path, ours: &Path) {
    let upstream = std::fs::read_to_string(upstream).unwrap();
    let ours = std::fs::read_to_string(ours).unwrap();
    let upstream = upstream
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect::<Vec<_>>();
    let ours = ours
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(ours.len(), upstream.len());
    for (ours, upstream) in ours.into_iter().zip(upstream) {
        let ours = ours.split('\t').collect::<Vec<_>>();
        let upstream = upstream.split('\t').collect::<Vec<_>>();
        assert_eq!(&ours[..3], &upstream[..3]);
        for (ours, upstream) in ours[3..].iter().zip(&upstream[3..]) {
            let ours = ours.parse::<f64>().unwrap();
            let upstream = upstream.parse::<f64>().unwrap();
            assert!((ours - upstream).abs() <= 1e-4, "{ours} != {upstream}");
        }
    }
}

#[test]
fn target_files_use_shared_coordinates_and_product_overlap_policy() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("signals.vcf");
    let targets = directory.path().join("targets.txt");
    write_fixture(&input);
    std::fs::write(&targets, "chr1\t14\n").unwrap();

    let selected = |overlap| {
        analyze_selected(
            &input,
            SampleSelection::default(),
            CallOptions::default(),
            &SiteSelection {
                targets_file: Some(targets.clone()),
                targets_overlap: overlap,
                ..SiteSelection::default()
            },
        )
        .map(|result| result.chromosomes[0].sites.len())
    };

    assert!(selected(OverlapMode::Position).is_err());
    assert_eq!(selected(OverlapMode::Record).unwrap(), 1);
    assert_eq!(selected(OverlapMode::Variant).unwrap(), 1);
}

#[test]
fn regions_require_an_index_and_compose_with_streaming_targets() {
    let directory = tempfile::tempdir().unwrap();
    let plain = directory.path().join("signals.vcf");
    let input = directory.path().join("signals.vcf.gz");
    let regions = directory.path().join("regions.bed");
    let targets = directory.path().join("targets.txt");
    write_fixture(&plain);
    write_indexed_vcf(&plain, &input);
    std::fs::write(&regions, "chr1\t0\t20\n").unwrap();
    std::fs::write(&targets, "chr1\t20\n").unwrap();
    let selection = SiteSelection {
        regions_file: Some(regions),
        targets_file: Some(targets),
        ..SiteSelection::default()
    };

    let error = analyze_selected(
        &plain,
        SampleSelection::default(),
        CallOptions::default(),
        &selection,
    )
    .unwrap_err();
    assert!(error.to_string().contains("indexed"), "{error}");

    let result = analyze_selected(
        &input,
        SampleSelection::default(),
        CallOptions::default(),
        &selection,
    )
    .unwrap();
    assert_eq!(result.chromosomes[0].sites.len(), 1);
    assert_eq!(result.chromosomes[0].sites[0].position, 20);
}

#[test]
fn indexed_regions_preserve_header_order_for_allele_frequency_merging() {
    let directory = tempfile::tempdir().unwrap();
    let plain = directory.path().join("signals.vcf");
    let input = directory.path().join("signals.vcf.gz");
    let regions = directory.path().join("regions.txt");
    let frequencies = directory.path().join("frequencies.tsv");
    write_fixture(&plain);
    write_indexed_vcf(&plain, &input);
    std::fs::write(&regions, "chr1\t20\nchr1\t10\n").unwrap();
    std::fs::write(
        &frequencies,
        "chr1\t10\tAACGT,AATGT\t0.2\nchr1\t20\tA,G\t0.3\n",
    )
    .unwrap();
    let result = analyze_with_allele_frequencies_selected(
        &input,
        SampleSelection::default(),
        CallOptions::default(),
        &frequencies,
        &SiteSelection {
            regions_file: Some(regions),
            ..SiteSelection::default()
        },
    )
    .unwrap();
    assert_eq!(
        result.chromosomes[0]
            .sites
            .iter()
            .map(|site| site.position)
            .collect::<Vec<_>>(),
        [10, 20]
    );
}

#[test]
fn inline_target_exclusion_uses_the_same_overlap_contract() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("signals.vcf");
    write_fixture(&input);
    let selection = SiteSelection {
        targets: vec!["chr1:10".to_owned()],
        exclude_targets: true,
        ..SiteSelection::default()
    };
    let result = analyze_selected(
        &input,
        SampleSelection::default(),
        CallOptions::default(),
        &selection,
    )
    .unwrap();
    assert_eq!(result.chromosomes[0].sites.len(), 1);
    assert_eq!(result.chromosomes[0].sites[0].position, 20);

    let error = analyze_selected(
        &input,
        SampleSelection::default(),
        CallOptions::default(),
        &SiteSelection {
            exclude_targets: true,
            ..SiteSelection::default()
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("requires a target"), "{error}");
}

#[test]
fn indexed_regions_suppress_records_spanning_multiple_queries() {
    let directory = tempfile::tempdir().unwrap();
    let plain = directory.path().join("spanning.vcf");
    let input = directory.path().join("spanning.vcf.gz");
    let regions = directory.path().join("regions.txt");
    std::fs::write(
        &plain,
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"BAF\">\n\
##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"LRR\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t10\t.\tA\t<DEL>\t.\t.\tEND=30\tBAF:LRR\t0.5:0.0\n",
    )
    .unwrap();
    write_indexed_vcf(&plain, &input);
    std::fs::write(&regions, "chr1\t11\nchr1\t20\n").unwrap();
    let result = analyze_selected(
        &input,
        SampleSelection::default(),
        CallOptions::default(),
        &SiteSelection {
            regions_file: Some(regions),
            ..SiteSelection::default()
        },
    )
    .unwrap();
    assert_eq!(result.chromosomes[0].sites.len(), 1);
}

#[test]
fn polysomy_applies_streaming_targets_before_histogramming() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("polysomy.vcf");
    let targets = directory.path().join("targets.txt");
    write_polysomy_fixture(&input);
    std::fs::write(&targets, "chr1\t1\t150\n").unwrap();
    let result = analyze_distributions_selected(
        &input,
        Some("SAMPLE".to_owned()),
        PolysomyOptions::default(),
        &SiteSelection {
            targets_file: Some(targets),
            ..SiteSelection::default()
        },
    )
    .unwrap();
    assert_eq!(result.chromosomes[0].observations, 150);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn target_overlap_modes_match_bcftools_1_24() {
    let version = Command::new("bcftools").arg("--version").output().unwrap();
    assert!(version.status.success());
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24"));
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("signals.vcf");
    let targets = directory.path().join("targets.txt");
    write_fixture(&input);
    std::fs::write(&targets, "chr1\t11\nchr1\t20\n").unwrap();

    for (value, overlap) in [
        ("0", OverlapMode::Position),
        ("1", OverlapMode::Record),
        ("2", OverlapMode::Variant),
    ] {
        let upstream_output = directory.path().join(format!("bcftools-{value}"));
        let upstream = Command::new("bcftools")
            .args(["cnv", "-s", "SAMPLE", "--targets-overlap", value, "-T"])
            .arg(&targets)
            .arg("-o")
            .arg(&upstream_output)
            .arg(&input)
            .output()
            .unwrap();
        assert!(
            upstream.status.success(),
            "{}",
            String::from_utf8_lossy(&upstream.stderr)
        );
        let result = analyze_selected(
            &input,
            SampleSelection::default(),
            CallOptions::default(),
            &SiteSelection {
                targets_file: Some(targets.clone()),
                targets_overlap: overlap,
                ..SiteSelection::default()
            },
        )
        .unwrap();
        let ours_output = directory.path().join(format!("rsomics-{value}"));
        write_call_reports(&ours_output, &result).unwrap();
        assert_eq!(
            std::fs::read(upstream_output.join("dat.SAMPLE.tab")).unwrap(),
            std::fs::read(ours_output.join("dat.SAMPLE.tab")).unwrap(),
            "overlap mode {value}: dat.SAMPLE.tab"
        );
        assert_cn_close(
            &upstream_output.join("cn.SAMPLE.tab"),
            &ours_output.join("cn.SAMPLE.tab"),
        );
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn indexed_region_overlap_modes_match_bcftools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let plain = directory.path().join("signals.vcf");
    let input = directory.path().join("signals.vcf.gz");
    let bcf = directory.path().join("signals.bcf");
    let regions = directory.path().join("regions.txt");
    write_fixture(&plain);
    write_indexed_vcf(&plain, &input);
    let conversion = Command::new("bcftools")
        .args(["view", "-Ob", "-o"])
        .arg(&bcf)
        .arg(&plain)
        .output()
        .unwrap();
    assert!(
        conversion.status.success(),
        "{}",
        String::from_utf8_lossy(&conversion.stderr)
    );
    let indexing = Command::new("bcftools")
        .args(["index", "--force"])
        .arg(&bcf)
        .output()
        .unwrap();
    assert!(
        indexing.status.success(),
        "{}",
        String::from_utf8_lossy(&indexing.stderr)
    );
    std::fs::write(&regions, "chr1\t11\nchr1\t20\n").unwrap();

    for (encoding, input) in [("vcf", &input), ("bcf", &bcf)] {
        for (value, overlap) in [
            ("0", OverlapMode::Position),
            ("1", OverlapMode::Record),
            ("2", OverlapMode::Variant),
        ] {
            let upstream_output = directory
                .path()
                .join(format!("bcftools-{encoding}-{value}"));
            let upstream = Command::new("bcftools")
                .args(["cnv", "-s", "SAMPLE", "--regions-overlap", value, "-R"])
                .arg(&regions)
                .arg("-o")
                .arg(&upstream_output)
                .arg(input)
                .output()
                .unwrap();
            assert!(
                upstream.status.success(),
                "{}",
                String::from_utf8_lossy(&upstream.stderr)
            );
            let result = analyze_selected(
                input,
                SampleSelection::default(),
                CallOptions::default(),
                &SiteSelection {
                    regions_file: Some(regions.clone()),
                    regions_overlap: overlap,
                    ..SiteSelection::default()
                },
            )
            .unwrap();
            let ours_output = directory.path().join(format!("rsomics-{encoding}-{value}"));
            write_call_reports(&ours_output, &result).unwrap();
            assert_eq!(
                std::fs::read(upstream_output.join("dat.SAMPLE.tab")).unwrap(),
                std::fs::read(ours_output.join("dat.SAMPLE.tab")).unwrap(),
                "{encoding} overlap mode {value}: dat.SAMPLE.tab"
            );
            assert_cn_close(
                &upstream_output.join("cn.SAMPLE.tab"),
                &ours_output.join("cn.SAMPLE.tab"),
            );
        }
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24 with polysomy"]
fn polysomy_region_and_target_selection_match_bcftools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let plain = directory.path().join("polysomy.vcf");
    let indexed = directory.path().join("polysomy.vcf.gz");
    let selection_file = directory.path().join("selection.txt");
    write_polysomy_fixture(&plain);
    write_indexed_vcf(&plain, &indexed);
    std::fs::write(&selection_file, "chr1\t1\t150\n").unwrap();

    for (kind, input) in [("targets", &plain), ("regions", &indexed)] {
        let upstream_output = directory.path().join(format!("bcftools-polysomy-{kind}"));
        let mut command = Command::new("bcftools");
        command.args(["polysomy", "-s", "SAMPLE"]);
        command.arg(if kind == "targets" { "-T" } else { "-R" });
        let upstream = command
            .arg(&selection_file)
            .arg("-o")
            .arg(&upstream_output)
            .arg(input)
            .output()
            .unwrap();
        assert!(
            upstream.status.success(),
            "{}",
            String::from_utf8_lossy(&upstream.stderr)
        );
        let selection = if kind == "targets" {
            SiteSelection {
                targets_file: Some(selection_file.clone()),
                ..SiteSelection::default()
            }
        } else {
            SiteSelection {
                regions_file: Some(selection_file.clone()),
                ..SiteSelection::default()
            }
        };
        let ours = analyze_polysomy_selected(
            input,
            Some("SAMPLE".to_owned()),
            PolysomyOptions::default(),
            &selection,
        )
        .unwrap();
        let ours_output = directory.path().join(format!("rsomics-polysomy-{kind}"));
        write_polysomy_reports(&ours_output, &ours).unwrap();
        let rows = |path: &Path, prefix: &str| {
            std::fs::read_to_string(path)
                .unwrap()
                .lines()
                .filter(|line| line.starts_with(prefix))
                .map(str::to_owned)
                .collect::<Vec<_>>()
        };
        assert_eq!(
            rows(&upstream_output.join("dist.dat"), "DIST\t"),
            rows(&ours_output.join("dist.dat"), "DIST\t"),
            "{kind} distribution"
        );
        let upstream_cn = rows(&upstream_output.join("dist.dat"), "CN\t");
        let ours_cn = rows(&ours_output.join("dist.dat"), "CN\t");
        assert_eq!(ours_cn.len(), upstream_cn.len());
        for (ours, upstream) in ours_cn.iter().zip(upstream_cn) {
            let ours = ours.split('\t').collect::<Vec<_>>();
            let upstream = upstream.split('\t').collect::<Vec<_>>();
            assert_eq!(&ours[..3], &upstream[..3], "{kind} copy number");
            let ours = ours[3].parse::<f64>().unwrap();
            let upstream = upstream[3].parse::<f64>().unwrap();
            assert!((ours - upstream).abs() <= 0.2, "{ours} != {upstream}");
        }
    }
}
