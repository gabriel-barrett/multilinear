use p3_air::{Air, AirBuilderWithPublicValues, BaseAir};
use p3_baby_bear;
use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_challenger::{DuplexChallenger, HashChallenger, SerializingChallenger32};
use p3_commit::ExtensionMmcs;
use p3_goldilocks::Goldilocks;
// use p3_dft::Radix2Dit as Radix2DitParallel;
use p3_dft::{Radix2DitParallel, TwoAdicSubgroupDft};
use p3_field::extension::BinomialExtensionField;
use p3_field::{Field, PrimeCharacteristicRing};
use p3_fri::{FriConfig, TwoAdicFriPcs};
use p3_keccak::{Keccak256Hash, KeccakF};
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::Matrix;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{
    CompressionFunctionFromHasher, PaddingFreeSponge, SerializingHasher, TruncatedPermutation,
};
use p3_uni_stark::{prove, verify, Proof, StarkConfig, StarkGenericConfig};
use rand::rngs::SmallRng;
use rand::SeedableRng;

use crate::benchmark;

pub const fn fri_config<Mmcs>(mmcs: Mmcs) -> FriConfig<Mmcs> {
    FriConfig {
        log_blowup: 1,
        log_final_poly_len: 0,
        num_queries: 128,
        proof_of_work_bits: 0,
        mmcs,
    }
}

struct CS {}

const LOG_WIDTH: usize = 6;
const WIDTH: usize = 1 << LOG_WIDTH;

impl<F> BaseAir<F> for CS {
    fn width(&self) -> usize {
        WIDTH
    }
}

impl<AB: AirBuilderWithPublicValues> Air<AB> for CS {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0).unwrap();
        assert!(LOG_WIDTH > 1);
        local.chunks_exact(4).for_each(|chunk| {
            builder.assert_eq(
                chunk[0] * chunk[0] + chunk[1] * chunk[1],
                chunk[2] * chunk[2],
            );
            builder.assert_eq(chunk[0] + chunk[1], chunk[3]);
        });
        // builder.assert_zero(AB::F::ZERO);
    }
}

fn demo_trace<F: Field>(total_log_height: u8) -> RowMajorMatrix<F> {
    let f = F::from_u32;
    let log_height = 10 - LOG_WIDTH as u8;
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
    for _ in 0..(4 + total_log_height - log_height) {
        trace.extend(trace.clone());
    }
    assert_eq!(trace.len(), WIDTH * (1 << total_log_height));
    RowMajorMatrix::new(trace, WIDTH)
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
type Dft<F> = Radix2DitParallel<F>;
type Pcs = TwoAdicFriPcs<Val, Dft<Val>, ValMmcs, ChallengeMmcs>;
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
const TOTAL_LOG_HEIGHT: u8 = 20;

const KECCAK_VECTOR_LEN: usize = p3_keccak::VECTOR_LEN;
type KeccakCompressionFunction =
    CompressionFunctionFromHasher<PaddingFreeSponge<KeccakF, 25, 17, 4>, 2, 4>;
type KeccakMerkleMmcs<F> = MerkleTreeMmcs<
    [F; KECCAK_VECTOR_LEN],
    [u64; KECCAK_VECTOR_LEN],
    SerializingHasher<PaddingFreeSponge<KeccakF, 25, 17, 4>>,
    KeccakCompressionFunction,
    4,
>;
pub(crate) type KeccakStarkConfig<F, EF, DFT> = StarkConfig<
    TwoAdicFriPcs<F, DFT, KeccakMerkleMmcs<F>, ExtensionMmcs<F, EF, KeccakMerkleMmcs<F>>>,
    EF,
    SerializingChallenger32<F, HashChallenger<u8, Keccak256Hash, 32>>,
>;

const fn get_keccak_mmcs<F: Field>() -> KeccakMerkleMmcs<F> {
    let u64_hash = PaddingFreeSponge::<KeccakF, 25, 17, 4>::new(KeccakF {});

    let field_hash = SerializingHasher::new(u64_hash);

    let compress = KeccakCompressionFunction::new(u64_hash);

    KeccakMerkleMmcs::new(field_hash, compress)
}

#[test]
fn p3_test() {
    let mut rng = SmallRng::seed_from_u64(1);
    let perm = Perm::new_from_rng_128(&mut rng);
    let hash = MyHash::new(perm.clone());
    let compress = MyCompress::new(perm.clone());
    let val_mmcs = ValMmcs::new(hash, compress);
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
    let dft = Dft::default();
    let trace = demo_trace(TOTAL_LOG_HEIGHT);
    let fri_config = fri_config(challenge_mmcs);
    let pcs = Pcs::new(dft, val_mmcs, fri_config);
    let challenger = Challenger::new(perm.clone());
    let config = MyConfig::new(pcs, challenger);
    let pis = vec![];
    println!(
        "GENERATING POLYNOMIAL FOR HEIGHT {} AND WIDTH {}",
        1 << TOTAL_LOG_HEIGHT,
        WIDTH,
    );
    let proof = benchmark!("Proof: ", prove(&config, &CS {}, trace, &pis));
    report_proof_size(&proof);
    benchmark!(
        "Verification ",
        verify(&config, &CS {}, &proof, &pis).expect("verification failed")
    );
}

#[test]
fn p3_ntt_test() {
    let dft = Dft::<BabyBear>::default();
    let trace = demo_trace(TOTAL_LOG_HEIGHT);
    println!(
        "DFT BABYBEAR FOR HEIGHT {} AND WIDTH {}",
        1 << TOTAL_LOG_HEIGHT,
        WIDTH,
    );
    benchmark!("DFT: ", dft.dft_batch(trace));
    let dft = Dft::<Goldilocks>::default();
    let trace = demo_trace(TOTAL_LOG_HEIGHT);
    println!(
        "DFT GOLDILOCKS FOR HEIGHT {} AND WIDTH {}",
        1 << TOTAL_LOG_HEIGHT,
        WIDTH,
    );
    benchmark!("DFT: ", dft.dft_batch(trace));
}

#[test]
fn p3_keccak_test() {
    let dft = Dft::default();

    let val_mmcs = get_keccak_mmcs();

    let challenge_mmcs = ExtensionMmcs::<Val, Challenge, _>::new(val_mmcs.clone());
    let fri_config = fri_config(challenge_mmcs);

    let trace = demo_trace(TOTAL_LOG_HEIGHT);

    let pcs = TwoAdicFriPcs::new(dft, val_mmcs, fri_config);
    let challenger = SerializingChallenger32::from_hasher(vec![], Keccak256Hash {});

    let config = KeccakStarkConfig::new(pcs, challenger);
    let pis = vec![];
    println!(
        "GENERATING POLYNOMIAL FOR HEIGHT {} AND WIDTH {}",
        1 << TOTAL_LOG_HEIGHT,
        WIDTH,
    );
    let proof = benchmark!("Proof: ", prove(&config, &CS {}, trace, &pis));
    report_proof_size(&proof);
    benchmark!(
        "Verification ",
        verify(&config, &CS {}, &proof, &pis).expect("verification failed")
    );
}
