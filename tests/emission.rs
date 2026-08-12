use rsomics_cnv::emission::{
    AlleleFrequencies, EmissionModel, EvidenceParameters, SampleParameters,
};
use rsomics_cnv::signals::Measurement;

#[test]
fn default_emissions_match_bcftools_1_24() {
    let model =
        EmissionModel::new(SampleParameters::default(), EvidenceParameters::default()).unwrap();
    let frequencies = AlleleFrequencies::default();
    let cases = [
        (
            Measurement {
                baf: Some(0.5),
                lrr: Some(0.0),
            },
            [
                0.0,
                0.0001,
                0.999_930_172_169_829_4,
                0.000_239_442_208_458_273_2,
            ],
        ),
        (
            Measurement {
                baf: Some(1.0 / 3.0),
                lrr: Some(0.3),
            },
            [
                0.0,
                0.000_100_000_000_015_786_41,
                0.000_378_837_062_638_049_3,
                0.999_760_402_012_810_4,
            ],
        ),
        (
            Measurement {
                baf: Some(0.05),
                lrr: Some(-0.45),
            },
            [
                0.0,
                0.353_291_489_361_555_8,
                0.259_232_815_636_053_44,
                0.258_823_454_780_921_67,
            ],
        ),
    ];

    for (measurement, expected) in cases {
        let actual = model.probabilities(measurement, frequencies).unwrap();
        for (actual, expected) in actual.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
        }
    }
}

#[test]
fn missing_signal_uses_the_bcftools_cn0_prior() {
    let model =
        EmissionModel::new(SampleParameters::default(), EvidenceParameters::default()).unwrap();
    assert_eq!(
        model
            .probabilities(Measurement::default(), AlleleFrequencies::default())
            .unwrap(),
        [0.5, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]
    );
}

#[test]
fn allele_frequency_and_parameter_boundaries_are_checked() {
    let frequencies = AlleleFrequencies::from_non_reference_frequency(0.2).unwrap();
    assert!((frequencies.rr - 0.64).abs() < 1e-15);
    assert!((frequencies.ra - 0.32).abs() < 1e-15);
    assert!((frequencies.aa - 0.04).abs() < 1e-15);
    assert!(AlleleFrequencies::from_non_reference_frequency(-0.1).is_err());
    assert!(AlleleFrequencies::new(0.0, 0.0, 0.0).is_err());
    assert!(AlleleFrequencies::new(0.5, f64::NAN, 0.5).is_err());

    assert!(
        EmissionModel::new(
            SampleParameters {
                baf_deviation: 0.0,
                ..SampleParameters::default()
            },
            EvidenceParameters::default(),
        )
        .is_err()
    );
    assert!(
        EmissionModel::new(
            SampleParameters::default(),
            EvidenceParameters {
                lrr_weight: 1.1,
                ..EvidenceParameters::default()
            },
        )
        .is_err()
    );
}

#[test]
fn lrr_is_optional_only_when_its_weight_is_zero() {
    let measurement = Measurement {
        baf: Some(0.5),
        lrr: None,
    };
    let frequencies = AlleleFrequencies::default();
    let required =
        EmissionModel::new(SampleParameters::default(), EvidenceParameters::default()).unwrap();
    assert!(required.probabilities(measurement, frequencies).is_err());

    let baf_only = EmissionModel::new(
        SampleParameters::default(),
        EvidenceParameters {
            lrr_weight: 0.0,
            ..EvidenceParameters::default()
        },
    )
    .unwrap();
    assert!(baf_only.probabilities(measurement, frequencies).is_ok());
}
