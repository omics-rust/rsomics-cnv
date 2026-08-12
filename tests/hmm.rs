use rsomics_cnv::hmm::{Hmm, phred_error_probability};

#[test]
fn inference_matches_the_bcftools_1_24_hmm_oracle() {
    let positions = [99, 199, 499, 500];
    let emissions = [
        [0.01, 0.20, 0.70, 0.09],
        [0.01, 0.70, 0.20, 0.09],
        [0.01, 0.80, 0.10, 0.09],
        [0.01, 0.10, 0.80, 0.09],
    ];
    let expected = [
        [
            0.002_904_407_141_440_417_4,
            0.135_054_490_435_530_03,
            0.835_443_360_769_083_7,
            0.026_597_741_653_945_96,
        ],
        [
            1.280_115_341_655_358_3e-05,
            0.420_798_296_876_149_56,
            0.574_844_830_028_829_8,
            0.004_344_071_941_604_094,
        ],
        [
            0.000_366_758_660_327_208,
            0.452_627_471_904_765_3,
            0.517_998_280_003_419_4,
            0.029_007_489_431_488_113,
        ],
        [
            0.000_388_315_833_681_164_55,
            0.448_678_398_407_692_97,
            0.521_750_641_005_775_8,
            0.029_182_644_752_849_93,
        ],
    ];

    let result = Hmm::single_sample(1e-3)
        .unwrap()
        .infer(&positions, &emissions)
        .unwrap();
    assert_eq!(result.path, [2, 2, 2, 2]);
    for (actual, expected) in result.posteriors().zip(expected) {
        for (actual, expected) in actual.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-14, "{actual} != {expected}");
        }
    }
}

#[test]
fn paired_prior_prefers_joint_diploid_state() {
    let emissions = [[1.0; 16]];
    let result = Hmm::paired_samples(1e-9, 0.5)
        .unwrap()
        .infer(&[10], &emissions)
        .unwrap();
    assert_eq!(result.path, [10]);
    let posterior = result.posterior(0).unwrap();
    assert_eq!(posterior.len(), 16);
    assert!(posterior[10] > posterior[2]);
    assert!(result.posterior(usize::MAX / result.states()).is_none());
}

#[test]
fn invalid_models_and_observations_fail() {
    assert!(Hmm::single_sample(-1.0).is_err());
    assert!(Hmm::single_sample(0.3).is_err());
    assert!(Hmm::paired_samples(1e-9, f64::NAN).is_err());
    assert!(Hmm::paired_samples(1e-9, 1.1).is_err());

    let model = Hmm::single_sample(1e-9).unwrap();
    let no_emissions: [[f64; 4]; 0] = [];
    assert!(model.infer(&[1], &no_emissions).is_err());
    assert!(model.infer(&[2, 1], &[[1.0; 4], [1.0; 4]]).is_err());
    assert!(model.infer(&[1], &[[0.0; 4]]).is_err());
    assert!(model.infer(&[1], &[[1.0, 1.0, f64::NAN, 1.0]]).is_err());
}

#[test]
fn quality_accepts_an_error_probability_once() {
    assert_eq!(phred_error_probability(0.0).unwrap(), 99.0);
    assert!((phred_error_probability(0.1).unwrap() - 10.0).abs() < 1e-3);
    assert_eq!(phred_error_probability(1.0).unwrap(), 0.0);
    assert!(phred_error_probability(-0.1).is_err());
    assert!(phred_error_probability(f64::NAN).is_err());
}
