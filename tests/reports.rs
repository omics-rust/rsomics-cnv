use std::fs;

use rsomics_cnv::call::{CallResult, ChromosomeCall, RegionCall, SiteCall};
use rsomics_cnv::polysomy::{
    CandidateFit, ChromosomeDistribution, ChromosomePolysomy, DistributionBin, FitCurve,
    PolysomyResult,
};
use rsomics_cnv::reports::{write_call_reports, write_polysomy_reports};
use rsomics_cnv::signals::Measurement;

fn result(control: bool) -> CallResult {
    CallResult {
        sample: "QUERY".to_owned(),
        control_sample: control.then(|| "CONTROL".to_owned()),
        chromosomes: vec![ChromosomeCall {
            reference_name: "chr1".to_owned(),
            sites: vec![SiteCall {
                position: 11,
                copy_number: 2,
                control_copy_number: control.then_some(1),
                state_probability: 0.91,
                posterior: [0.01, 0.02, 0.94, 0.03],
                control_posterior: control.then_some([0.01, 0.91, 0.06, 0.02]),
                measurement: Measurement {
                    baf: Some(0.52),
                    lrr: Some(-0.12),
                },
                modeled_lrr: Some(-0.08),
                control_measurement: control.then_some(Measurement {
                    baf: Some(0.48),
                    lrr: Some(-0.31),
                }),
                control_modeled_lrr: control.then_some(-0.27),
            }],
            regions: vec![RegionCall {
                start: 11,
                end: 11,
                copy_number: 2,
                control_copy_number: control.then_some(1),
                quality: 10.4,
                sites: 1,
                heterozygous_sites: 1,
                control_heterozygous_sites: control.then_some(1),
            }],
        }],
    }
}

fn polysomy_result() -> PolysomyResult {
    PolysomyResult {
        sample: "SAMPLE".to_owned(),
        chromosomes: vec![ChromosomePolysomy {
            distribution: ChromosomeDistribution {
                reference_name: "chr1".to_owned(),
                observations: 200,
                bins: (0..150)
                    .map(|index| DistributionBin {
                        baf: index as f64 / 149.0,
                        normalized_count: if index == 74 { 1.0 } else { 0.0 },
                    })
                    .collect(),
                rr_boundary: 3,
                heterozygous_center: 75,
                aa_boundary: 146,
                fitted_end: 149,
                preliminary_copy_number: None,
            },
            copy_number: 2.0,
            absolute_deviation: Some(0.125),
            candidates: vec![CandidateFit {
                model_copy_number: 2,
                estimated_copy_number: 2.0,
                absolute_deviation: Some(0.125),
                selected: true,
                rejection: None,
                curves: vec![FitCurve {
                    start_bin: 3,
                    end_bin: 146,
                    absolute_deviation: 0.125,
                    function: "1.000000**2 * exp(-(x-0.500000)**2/0.040000**2)".to_owned(),
                }],
            }],
        }],
    }
}

#[test]
fn single_sample_report_bundle_is_complete() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("reports");
    write_call_reports(&output, &result(false)).unwrap();

    let mut names = fs::read_dir(&output)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect::<Vec<_>>();
    names.sort();
    assert_eq!(
        names,
        [
            "cn.QUERY.tab",
            "dat.QUERY.tab",
            "result.json",
            "summary.QUERY.tab"
        ]
    );

    let dat = fs::read_to_string(output.join("dat.QUERY.tab")).unwrap();
    assert!(dat.contains("chr1\t11\t0.520\t-0.120\n"));
    let cn = fs::read_to_string(output.join("cn.QUERY.tab")).unwrap();
    assert!(cn.contains("chr1\t11\t2\t0.010000\t0.020000\t0.940000\t0.030000\n"));
    let summary = fs::read_to_string(output.join("summary.QUERY.tab")).unwrap();
    assert!(summary.contains("RG\tchr1\t11\t11\t2\t10.4\t1\t1\n"));
    let json: serde_json::Value =
        serde_json::from_slice(&fs::read(output.join("result.json")).unwrap()).unwrap();
    assert_eq!(json["schema"], "rsomics-cnv/call-result/v1");
    assert_eq!(
        json["result"]["chromosomes"][0]["sites"][0]["modeled_lrr"],
        -0.08
    );
}

#[test]
fn paired_bundle_has_per_sample_and_joint_reports() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("reports");
    write_call_reports(&output, &result(true)).unwrap();

    let mut names = fs::read_dir(&output)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect::<Vec<_>>();
    names.sort();
    assert_eq!(
        names,
        [
            "cn.CONTROL.tab",
            "cn.QUERY.tab",
            "dat.CONTROL.tab",
            "dat.QUERY.tab",
            "result.json",
            "summary.CONTROL.tab",
            "summary.QUERY.tab",
            "summary.tab"
        ]
    );
    let joint = fs::read_to_string(output.join("summary.tab")).unwrap();
    assert!(joint.contains("RG\tchr1\t11\t11\t2\t1\t10.4\t1\t1\t1\t1\n"));
}

#[test]
fn failed_bundle_leaves_no_destination() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("reports");
    let mut invalid = result(false);
    invalid.sample = "../escape".to_owned();
    let error = write_call_reports(&output, &invalid).unwrap_err();
    assert!(error.to_string().contains("sample name"), "{error}");
    assert!(!output.exists());
    assert!(!directory.path().join("escape").exists());
}

#[test]
fn existing_destination_is_not_replaced() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("reports");
    fs::create_dir(&output).unwrap();
    fs::write(output.join("keep"), "original").unwrap();
    let error = write_call_reports(&output, &result(false)).unwrap_err();
    assert!(error.to_string().contains("already exists"), "{error}");
    assert_eq!(fs::read_to_string(output.join("keep")).unwrap(), "original");
}

#[test]
fn polysomy_report_bundle_contains_compatibility_and_machine_results() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("reports");
    write_polysomy_reports(&output, &polysomy_result()).unwrap();
    let mut names = fs::read_dir(&output)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect::<Vec<_>>();
    names.sort();
    assert_eq!(names, ["dist.dat", "result.json"]);
    let report = fs::read_to_string(output.join("dist.dat")).unwrap();
    assert_eq!(
        report
            .lines()
            .filter(|line| line.starts_with("DIST\t"))
            .count(),
        150
    );
    assert!(report.contains(
        "FIT\tchr1\t1.250000e-01\t3\t146\t1.000000**2 * exp(-(x-0.500000)**2/0.040000**2)"
    ));
    assert!(report.contains("CN\tchr1\t2.00\t0.125000\n"));
    let json: serde_json::Value =
        serde_json::from_slice(&fs::read(output.join("result.json")).unwrap()).unwrap();
    assert_eq!(json["schema"], "rsomics-cnv/polysomy-result/v1");
    assert_eq!(json["result"]["chromosomes"][0]["copy_number"], 2.0);
}

#[test]
fn invalid_polysomy_result_does_not_create_destination() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("reports");
    let mut result = polysomy_result();
    result.chromosomes[0].copy_number = f64::NAN;
    let error = write_polysomy_reports(&output, &result).unwrap_err();
    assert!(error.to_string().contains("copy number"), "{error}");
    assert!(!output.exists());
}
