use std::time::Instant;

use p3_air::{Air, AirBuilderWithPublicValues, BaseAir};
use p3_baby_bear;
use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{Field, PrimeCharacteristicRing};
use p3_fri::{create_benchmark_fri_config, TwoAdicFriPcs};
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::Matrix;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::{prove, verify, Proof, StarkConfig, StarkGenericConfig};
use rand::rngs::SmallRng;
use rand::SeedableRng;

struct PythagoreanCS {}

const PYTHAGOREAN_WIDTH: usize = 4;

impl<F> BaseAir<F> for PythagoreanCS {
    fn width(&self) -> usize {
        PYTHAGOREAN_WIDTH
    }
}

impl<AB: AirBuilderWithPublicValues> Air<AB> for PythagoreanCS {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0);
        let a = local[0];
        let b = local[1];
        let c = local[2];
        let d = local[3];

        builder.assert_zero(a * a + b * b - c * c);
        builder.assert_zero((a + b) * (a + b) - d * d);
    }
}

fn demo_trace(total_log_height: u8) -> RowMajorMatrix<BabyBear> {
    let f = BabyBear::from_u32;
    let log_height = 4;
    let mut trace = [
        3, 4, 5, 7, //
        5, 12, 13, 17, //
        8, 15, 17, 23, //
        7, 24, 25, 31, //
        20, 21, 29, 41, //
        12, 35, 37, 47, //
        9, 40, 41, 49, //
        28, 45, 53, 73, //
        11, 60, 61, 71, //
        16, 63, 65, 79, //
        33, 56, 65, 89, //
        48, 55, 73, 103, //
        13, 84, 85, 97, //
        36, 77, 85, 113, //
        39, 80, 89, 119, //
        65, 72, 97, 137, //
    ]
    .map(f)
    .to_vec();
    for _ in 0..(total_log_height - log_height) {
        trace.extend(trace.clone());
    }
    assert_eq!(trace.len(), PYTHAGOREAN_WIDTH * (1 << total_log_height));
    RowMajorMatrix::new(trace, PYTHAGOREAN_WIDTH)
}

type Val = BabyBear;
type Perm = Poseidon2BabyBear<16>;
type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
type ValMmcs =
    MerkleTreeMmcs<<Val as Field>::Packing, <Val as Field>::Packing, MyHash, MyCompress, 8>;
type Challenge = BinomialExtensionField<Val, 4>;
type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
type Challenger = DuplexChallenger<Val, Perm, 16, 8>;
type Dft = Radix2DitParallel<Val>;
type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;

#[inline]
pub fn report_proof_size<SC>(proof: &Proof<SC>)
where
    SC: StarkGenericConfig,
{
    let config = bincode::config::standard()
        .with_little_endian()
        .with_fixed_int_encoding();
    let proof_bytes =
        bincode::serde::encode_to_vec(proof, config).expect("Failed to serialize proof");
    println!("Proof size: {} bytes", proof_bytes.len());
}

#[test]
fn test_public_value_impl() {
    const TOTAL_LOG_HEIGHT: u8 = 20;
    let mut rng = SmallRng::seed_from_u64(1);
    let perm = Perm::new_from_rng_128(&mut rng);
    let hash = MyHash::new(perm.clone());
    let compress = MyCompress::new(perm.clone());
    let val_mmcs = ValMmcs::new(hash, compress);
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
    let dft = Dft::default();
    let trace = demo_trace(TOTAL_LOG_HEIGHT);
    let fri_config = create_benchmark_fri_config(challenge_mmcs);
    let pcs = Pcs::new(dft, val_mmcs, fri_config);
    let config = MyConfig::new(pcs);
    let mut challenger = Challenger::new(perm.clone());
    let pis = vec![];
    println!(
        "GENERATING POLYNOMIAL FOR HEIGHT {} AND WIDTH {}",
        1 << TOTAL_LOG_HEIGHT,
        PYTHAGOREAN_WIDTH,
    );
    let now = Instant::now();
    let proof = prove(&config, &PythagoreanCS {}, &mut challenger, trace, &pis);

    println!("Proof took {:?}", now.elapsed());
    report_proof_size(&proof);
    println!("VERIFYING");
    let now = Instant::now();
    let mut challenger = Challenger::new(perm);
    verify(&config, &PythagoreanCS {}, &mut challenger, &proof, &pis).expect("verification failed");
    println!("Verification took {:?}", now.elapsed());
}
