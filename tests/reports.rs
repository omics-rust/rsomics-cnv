use std::fs;

use rsomics_cnv::call::{AberrantEstimate, CallResult, ChromosomeCall, RegionCall, SiteCall};
use rsomics_cnv::polysomy::{
    CandidateFit, ChromosomeDistribution, ChromosomePolysomy, DistributionBin, FitCurve,
    FitRejection, PolysomyResult,
};
use rsomics_cnv::reports::{
    CallReportOptions, PolysomyReportOptions, write_call_reports, write_call_reports_with_options,
    write_polysomy_reports, write_polysomy_reports_with_options,
};
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
            query_estimate: None,
            control_estimate: None,
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
                    fitted: (3..=146)
                        .map(|index| DistributionBin {
                            baf: index as f64 / 149.0,
                            normalized_count: if index == 74 { 1.0 } else { 0.0 },
                        })
                        .collect(),
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
    assert_eq!(json["schema"], "rsomics-cnv/call-result/v2");
    assert_eq!(json["result"]["chromosomes"][0]["sites"], 1);
    assert_eq!(
        json["result"]["chromosomes"][0]["regions"][0]["copy_number"],
        2
    );
    assert_eq!(
        json["result"]["artifacts"]["query"]["data"],
        "dat.QUERY.tab"
    );
    assert!(json["result"]["chromosomes"][0]["posterior"].is_null());
}

#[test]
fn optimized_call_reports_include_checked_estimates() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("optimized");
    let mut result = result(false);
    result.chromosomes[0].query_estimate = Some(AberrantEstimate {
        fraction: 0.5,
        baf_deviation: 0.028_284,
    });
    write_call_reports(&output, &result).unwrap();
    let summary = std::fs::read_to_string(output.join("summary.QUERY.tab")).unwrap();
    assert!(summary.contains("CF\tchr1\t11\t11\t0.50\t0.028284\n"));
    let json: serde_json::Value =
        serde_json::from_slice(&std::fs::read(output.join("result.json")).unwrap()).unwrap();
    assert_eq!(
        json["result"]["chromosomes"][0]["query_estimate"]["fraction"],
        0.5
    );

    let invalid_output = directory.path().join("invalid-estimate");
    result.chromosomes[0].query_estimate = Some(AberrantEstimate {
        fraction: f64::NAN,
        baf_deviation: 0.028_284,
    });
    let error = write_call_reports(&invalid_output, &result).unwrap_err();
    assert!(error.to_string().contains("estimate"), "{error}");
    assert!(!invalid_output.exists());
}

#[test]
fn call_plot_is_a_transactional_self_contained_svg() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("call-plot");
    write_call_reports_with_options(
        &output,
        &result(false),
        CallReportOptions {
            plot_threshold: Some(0.0),
        },
    )
    .unwrap();
    let plot = std::fs::read_to_string(output.join("plot.QUERY.chr1.svg")).unwrap();
    assert!(plot.starts_with("<svg"), "{plot}");
    for label in ["LRR", "BAF", "Copy number", "chr1"] {
        assert!(plot.contains(label), "missing {label}: {plot}");
    }

    let invalid = directory.path().join("invalid-plot");
    let error = write_call_reports_with_options(
        &invalid,
        &result(false),
        CallReportOptions {
            plot_threshold: Some(f64::NAN),
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("plot threshold"), "{error}");
    assert!(!invalid.exists());
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

    let empty_output = directory.path().join("empty-regions");
    let mut empty_regions = result(false);
    empty_regions.chromosomes[0].regions.clear();
    let error = write_call_reports(&empty_output, &empty_regions).unwrap_err();
    assert!(error.to_string().contains("sites or regions"), "{error}");
    assert!(!empty_output.exists());
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
    assert!(json["result"]["chromosomes"][0]["candidates"][0]["curves"][0]["fitted"].is_null());
}

#[test]
fn polysomy_plots_include_distributions_and_copy_numbers() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("polysomy-plots");
    write_polysomy_reports_with_options(
        &output,
        &polysomy_result(),
        PolysomyReportOptions { plots: true },
    )
    .unwrap();
    let distribution = std::fs::read_to_string(output.join("distribution.chr1.svg")).unwrap();
    assert!(distribution.starts_with("<svg"), "{distribution}");
    assert!(distribution.contains("BAF distribution"), "{distribution}");
    let copy_number = std::fs::read_to_string(output.join("copy-number.svg")).unwrap();
    assert!(copy_number.starts_with("<svg"), "{copy_number}");
    assert!(copy_number.contains("Copy number"), "{copy_number}");

    let unresolved_output = directory.path().join("unresolved-plots");
    let mut unresolved = polysomy_result();
    unresolved.chromosomes[0].copy_number = -1.0;
    unresolved.chromosomes[0].absolute_deviation = None;
    unresolved.chromosomes[0].candidates[0].selected = false;
    unresolved.chromosomes[0].candidates[0].rejection = Some(FitRejection::FitThreshold);
    write_polysomy_reports_with_options(
        &unresolved_output,
        &unresolved,
        PolysomyReportOptions { plots: true },
    )
    .unwrap();
    let distribution =
        std::fs::read_to_string(unresolved_output.join("distribution.chr1.svg")).unwrap();
    assert!(distribution.contains("unresolved"), "{distribution}");
    let copy_number = std::fs::read_to_string(unresolved_output.join("copy-number.svg")).unwrap();
    assert!(copy_number.contains("Unresolved"), "{copy_number}");
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
