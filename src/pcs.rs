/*
use std::time::Instant;

use spongefish::DomainSeparator;
use spongefish_pow::blake3::Blake3PoW;
use whir::whir::{
    committer::CommitmentWriter, domainsep::WhirDomainSeparator, parameters::WhirConfig,
    prover::Prover, statement::Statement, verifier::Verifier, whir_proof_size,
};
use whir::{
    crypto::merkle_tree::{keccak::KeccakMerkleTreeParams, HashCounter},
    parameters::{
        default_max_pow, FoldType, FoldingFactor, MultivariateParameters, SoundnessType,
        WhirParameters,
    },
    poly_utils::{coeffs::CoefficientList, multilinear::MultilinearPoint},
    whir::{
        committer::CommitmentReader,
        statement::{StatementVerifier, Weights},
    },
};

use crate::{BaseField, ExtField};

pub type PowStrategy = Blake3PoW;
pub type MerkleConfig = KeccakMerkleTreeParams<ExtField>;

#[derive(Clone, Copy)]
pub struct PCSConfig {
    security_level: usize,
    rate: usize,
    first_round_folding_factor: usize,
    folding_factor: usize,
    soundness_type: SoundnessType,
    fold_optimisation: FoldType,
}

pub fn run_pcs(
    conf: PCSConfig,
    polynomial: CoefficientList<BaseField>,
    points: &[MultilinearPoint<ExtField>],
) {
    let security_level = conf.security_level;
    let starting_rate = conf.rate;
    let first_round_folding_factor = conf.first_round_folding_factor;
    let folding_factor = conf.folding_factor;
    let fold_optimisation = conf.fold_optimisation;
    let soundness_type = conf.soundness_type;

    let num_variables = polynomial.num_variables();
    let pow_bits = default_max_pow(num_variables, starting_rate);

    let mv_params = MultivariateParameters::<ExtField>::new(num_variables);
    let whir_params = WhirParameters::<MerkleConfig, PowStrategy> {
        initial_statement: true,
        security_level,
        pow_bits,
        folding_factor: FoldingFactor::ConstantFromSecondRound(
            first_round_folding_factor,
            folding_factor,
        ),
        leaf_hash_params: (),
        two_to_one_params: (),
        soundness_type,
        fold_optimisation,
        _pow_parameters: Default::default(),
        starting_log_inv_rate: starting_rate,
    };
    let params = WhirConfig::<ExtField, MerkleConfig, PowStrategy>::new(mv_params, whir_params);

    let domainsep = DomainSeparator::new("🌪️")
        .commit_statement(&params)
        .add_whir_proof(&params);

    let mut prover_state = domainsep.to_prover_state();

    if !params.check_pow_bits() {
        println!("WARN: more PoW bits required than what specified.");
    }

    let whir_execution_time = Instant::now();
    let mut statement: Statement<ExtField> = Statement::<ExtField>::new(num_variables);
    for point in points {
        let eval = polynomial.evaluate_at_extension(point);
        let weights = Weights::evaluation(point.clone());
        statement.add_constraint(weights, eval);
    }
    println!("Execution time: {:.1?}", whir_execution_time.elapsed());

    let whir_commitment_time = Instant::now();
    let committer = CommitmentWriter::new(params.clone());
    let witness = committer.commit(&mut prover_state, polynomial).unwrap();
    println!("Commitment time: {:.1?}", whir_commitment_time.elapsed());

    let whir_prover_time = Instant::now();
    let prover = Prover(params.clone());
    let proof = prover
        .prove(&mut prover_state, statement.clone(), witness)
        .unwrap();
    println!("Prover time: {:.1?}", whir_prover_time.elapsed());
    println!(
        "Proof size: {:.1} KiB",
        whir_proof_size(prover_state.narg_string(), &proof) as f64 / 1024.0
    );

    let statement_verifier = StatementVerifier::from_statement(&statement);
    // Just not to count that initial inversion (which could be precomputed)
    let commitment_reader = CommitmentReader::new(&params);
    let verifier = Verifier::new(&params);
    HashCounter::reset();
    let whir_verifier_time = Instant::now();
    let mut verifier_state = domainsep.to_verifier_state(prover_state.narg_string());
    let parsed_commitment = commitment_reader
        .parse_commitment(&mut verifier_state)
        .unwrap();
    verifier
        .verify(
            &mut verifier_state,
            &parsed_commitment,
            &statement_verifier,
            &proof,
        )
        .unwrap();
    println!("Verifier time: {:.1?}", whir_verifier_time.elapsed());
    println!(
        "Average hashes: {:.1}k",
        (HashCounter::get() as f64) / 1000.0
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_conf() -> PCSConfig {
        PCSConfig {
            security_level: 100,
            rate: 1,
            first_round_folding_factor: 4,
            folding_factor: 4,
            soundness_type: SoundnessType::UniqueDecoding,
            fold_optimisation: FoldType::ProverHelps,
        }
    }

    fn default_polynomial(num_variables: usize) -> CoefficientList<BaseField> {
        let num_coeffs = 1 << num_variables;
        CoefficientList::new((0..num_coeffs).map(BaseField::from).collect())
    }

    fn default_points(
        num_variables: usize,
        num_evaluations: usize,
    ) -> Vec<MultilinearPoint<ExtField>> {
        (0..num_evaluations)
            .map(|x| MultilinearPoint(vec![ExtField::from(x as u64); num_variables]))
            .collect()
    }

    #[test]
    fn pcs_test() {
        let num_variables = 10;
        let polynomial = default_polynomial(num_variables);
        let conf = default_conf();
        let num_evaluations = 1;
        let points = default_points(num_variables, num_evaluations);
        run_pcs(conf, polynomial, &points);
    }
}
*/
