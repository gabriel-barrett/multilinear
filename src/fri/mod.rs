use crate::field::Field128;
use crate::merkle_tree::{Merkle, MerkleInclusionPath, MerkleInclusionPathError};
use crate::ntt::{NttField, Polynomial};
use crate::{field::Field, merkle_tree::HashDigest};
use bincode; // Ensure bincode is imported
use serde::{Deserialize, Serialize};

use sha2::{Digest, Sha256};

pub struct Transcript {
    state: Sha256,
}

impl Default for Transcript {
    fn default() -> Self {
        Self::new()
    }
}

impl Transcript {
    pub fn new() -> Self {
        Self {
            state: Sha256::new(),
        }
    }

    pub fn append_message(&mut self, label: &[u8], message: &[u8]) {
        self.state.update(label);
        self.state.update(message);
    }

    pub fn random(&self) -> [u8; 32] {
        let cloned_state = self.state.clone();
        let result = cloned_state.finalize();
        let mut random_bytes = [0u8; 32];
        random_bytes.copy_from_slice(&result[..32]);
        random_bytes
    }
}

pub trait HashableField: Field + AsRef<[u8]> {
    fn to_bytes(&self) -> &[u8];
    fn from_digest(digest: &[u8; 32]) -> Self;
}

impl HashableField for Field128 {
    fn to_bytes(&self) -> &[u8] {
        self.as_ref()
    }

    fn from_digest(digest: &[u8; 32]) -> Self {
        let x = u128::from_le_bytes(digest[0..16].try_into().unwrap());
        Self::from(x)
    }
}

pub struct ProverData<F> {
    pub commitments: Vec<Merkle<F>>,
    pub polynomials: Vec<Vec<F>>,
    pub last_element: Option<F>,
}

pub const LOG_BLOWUP: usize = 1;
pub const NUM_QUERIES: usize = 128;

pub fn reed_solomon<F: Field>(mut coeffs: Vec<F>, gen: F) -> Vec<F> {
    // first, multiply the size of `coeffs` by a factor of `blowup` through adding zeros
    let n = coeffs.len();
    let blowup = 1 << LOG_BLOWUP;
    assert!(blowup > 1);
    coeffs.resize(blowup * n, F::from(0));
    // use `ntt` to compute the Reed-Solomon encoding.
    let lagrange = Polynomial { coeffs }.ntt_iterative(gen);
    lagrange.evals
}

impl<F: HashableField> ProverData<F> {
    pub fn init(values: Vec<F>, gen: F, transcript: &mut Transcript) -> Self {
        // `values` must be power of two.
        assert!(
            values.len().is_power_of_two(),
            "Input size must be a power of two"
        );
        // push save a copy of `values` to `polynomials`
        let polynomials = vec![values.clone()];
        // use `reed_solomon` to compute the values for commitment.
        let rs_encoded = reed_solomon(values, gen);
        // commit to a `Merkle` tree using `to_bytes` method.
        let mut commitments = Vec::new();
        let merkle = Merkle::commit(rs_encoded);
        let root = merkle.root();
        // add to `commitments`.
        commitments.push(merkle);
        // Use the `root()` to update the transcript
        transcript.append_message(b"merkle_root", root.as_slice());
        Self {
            commitments,
            polynomials,
            last_element: None,
        }
    }

    pub fn fold_step(&mut self, gen: F, transcript: &mut Transcript) {
        let last_poly = self.polynomials.last().unwrap().clone();
        let n = last_poly.len();
        if n <= 1 {
            return;
        }

        // generate random field element called `r` from the transcript using `random` and `from_digest`
        let random_bytes = transcript.random();
        let r = F::from_digest(&random_bytes);

        let half_n = n >> 1;
        let mut next_poly = Vec::with_capacity(half_n);

        for i in 0..(half_n) {
            let even = last_poly[i * 2];
            let odd = last_poly[i * 2 + 1];

            next_poly.push(even + r * odd);
        }
        if half_n == 1 {
            // sanity check: last polynomial must be constant
            let first = next_poly[0];
            assert!(
                next_poly.iter().all(|next| first == *next),
                "not an RS code"
            );
            self.last_element = Some(first);
            transcript.append_message(b"last_element", first.as_ref());
            return;
        }
        self.polynomials.push(next_poly.clone());

        // Use `reed_solomon` to compute the values for commitment.
        let next_gen = gen * gen;
        let rs_encoded = reed_solomon(next_poly, next_gen);

        // `commit` to Merkle, etc
        let merkle = Merkle::commit(rs_encoded);
        let root = merkle.root();
        self.commitments.push(merkle);

        // Use the `root()` to update the transcript
        transcript.append_message(b"merkle_root", root.as_slice());
    }

    pub fn fold_step_opt(&mut self, gen: F, transcript: &mut Transcript) {
        // do not use polynomials. instead work solely on lagrange basis reading
        // the merkle `value` field
        let last_data = self.commitments.last().unwrap().data.clone();
        let n = last_data.len();
        let blowup = 1 << LOG_BLOWUP;
        if n <= blowup {
            return;
        }
        let random_bytes = transcript.random();
        let r = F::from_digest(&random_bytes);
        let half_n = n >> 1;
        let mut next_data = Vec::with_capacity(half_n);
        let mut gen_pow = F::from(1);
        for i in 0..half_n {
            // p(gen^i)
            let a = last_data[i];
            // p(-gen^i)
            let b = last_data[i + half_n];
            // even(x^2) = (p(x) + p(-x))/2, where x = gen^i
            let even = (a + b) / F::from(2);
            // odd(x^2) = (p(x) - p(-x))/2x, where x = gen^i
            let odd = (a - b) / (F::from(2) * gen_pow);
            // p(x) + p(-x) == 2*even(x^2)
            next_data.push(even + r * odd);
            gen_pow *= gen;
        }

        if half_n == blowup {
            // sanity check: last RS code must be constant
            let first = next_data[0];
            assert!(
                next_data.iter().all(|next| first == *next),
                "not an RS code"
            );
            self.last_element = Some(first);
            transcript.append_message(b"last_element", first.as_ref());
            return;
        }
        // `commit` to Merkle, etc
        let merkle = Merkle::commit(next_data);
        let root = merkle.root();
        self.commitments.push(merkle);

        // Use the `root()` to update the transcript
        transcript.append_message(b"merkle_root", root.as_slice());
    }

    pub fn fold(gen: F, values: Vec<F>, transcript: &mut Transcript) -> Self {
        let mut prover_data = Self::init(values, gen, transcript);
        let mut current_gen = gen;
        while prover_data.last_element.is_none() {
            prover_data.fold_step_opt(current_gen, transcript);
            current_gen *= current_gen;
        }
        prover_data
    }

    pub fn fold_roots(&self) -> Vec<HashDigest> {
        self.commitments
            .iter()
            .map(|merkle| merkle.root())
            .collect()
    }

    pub fn open_query_at(&self, index: usize) -> QueryProof<F> {
        let n = self.commitments[0].data.len();
        assert!(index < n / 2);
        let conjugate_index = index + n / 2;

        let mut paths = Vec::new();
        let mut current_index = index;
        let mut current_conjugate = conjugate_index;

        for merkle in &self.commitments {
            let path = merkle.open(current_index).expect("Index out of bounds");
            let conjugate_path = merkle
                .open(current_conjugate)
                .expect("Conjugate index out of bounds");

            paths.push((path, conjugate_path));

            current_index /= 2;
            current_conjugate /= 2;
        }

        QueryProof { index, paths }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryProof<F> {
    // initial random index, from 0..N/2
    pub index: usize,
    // merkle paths for all fold layers at both the index and index + N/2
    // the index at subsequent layers are halved
    pub paths: Vec<(MerkleInclusionPath<F>, MerkleInclusionPath<F>)>,
}

impl<F: HashableField + NttField> QueryProof<F> {
    pub fn verify(
        &self,
        commitments: &[HashDigest],
        gen: F,
        last_element: F,
        random_elements: &[F],
    ) -> Result<(), FriProofError> {
        if self.paths.len() != commitments.len() {
            return Err(FriProofError::WrongNumberOfPaths);
        }
        for ((value_path, minus_value_path), commitment) in
            self.paths.iter().zip(commitments.iter())
        {
            if let Err(err) = value_path.verify(commitment) {
                return Err(FriProofError::InclusionPathError(err));
            }
            if let Err(err) = minus_value_path.verify(commitment) {
                return Err(FriProofError::InclusionPathError(err));
            }
        let mut current_gen = gen;
        for i in 0..self.paths.len() {
            let (value_path, minus_value_path) = &self.paths[i];
            let value = value_path.value; // p(g^i)
            let minus_value = minus_value_path.value; // p(-g^i)

            let even = (value + minus_value) / F::from(2);
            let odd = (value - minus_value) / (F::from(2) * current_gen);

            if i == self.paths.len() - 1 {
                // Último caso: verificar com o elemento final
                if last_element != even + random_elements[i] * odd {
                    return Err(FriProofError::QueryMismatch(i));
                }
            } else {
                let (next_value_path, _) = &self.paths[i + 1];
                let next_value = next_value_path.value; // p(g^(2i))
                if next_value != even + random_elements[i] * odd {
                    return Err(FriProofError::QueryMismatch(i));
                }
                current_gen *= current_gen;
            }
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FriProof<F: HashableField + NttField> {
    // there are N commitments, where N is the log2 of the message size (the polynomial)
    pub commitments: Vec<HashDigest>,
    // there are `NUM_QUERIES` number of query proofs
    pub queries: Vec<QueryProof<F>>,
    // this is the last message, a single element
    pub last_elem: F,
}

pub enum FriProofError {
    QueryMismatch(usize),
    WrongNumberOfQueries,
    WrongNumberOfPaths,
    InclusionPathError(MerkleInclusionPathError),
    Generic,
}

impl<F: HashableField + NttField> FriProof<F> {
    pub fn prove(message: Vec<F>, transcript: &mut Transcript) -> FriProof<F> {
        // get the generator for length = blowup * message.len
        let n = message.len();
        let blowup = 1 << LOG_BLOWUP;
        let domain_size = blowup * n;
        let log_size = domain_size.trailing_zeros();
        let gen = F::pow_2_generator(log_size as u64).unwrap();

        // call `fold`
        let prover_data = ProverData::fold(gen, message, transcript);
        // for `0..NUM_QUERIES` generate random index between `0..domain_size/2`
        let mut queries = Vec::with_capacity(NUM_QUERIES);
        for _ in 0..NUM_QUERIES {
            let random_bytes = transcript.random();
            let random_index = (u64::from_le_bytes(random_bytes[..8].try_into().unwrap())
                % (domain_size / 2) as u64) as usize;
            // open query at this index and add the proof to a vector of query proofs
            let query_proof = prover_data.open_query_at(random_index);
            queries.push(query_proof);
            // use the `index` to update the transcript
            transcript.append_message(b"query_index", &random_index.to_le_bytes());
        }
        // at the end create the FriProof using the queries, last_elem and the
        FriProof {
            commitments: prover_data.fold_roots(),
            queries,
            last_elem: prover_data.last_element.unwrap(),
        }
    }

    pub fn verify(&self) -> Result<(), FriProofError> {
        // `verify` has to simulate two stages, namely the "fold" stage and the "query" stage.  The
        // "fold" stage will basically just produce a bunch of random values to pass to the query
        // verifier The "query" stage will call the query verifier for all queries of the proof Also
        // verify that the number of queries is equal to `NUM_QUERIES`
        if self.queries.len() != NUM_QUERIES {
            return Err(FriProofError::WrongNumberOfQueries);
        }

        // Create a transcript for verification
        let mut transcript = Transcript::new();
        let mut random_elements = Vec::new();
        let mut current_gen =
            F::pow_2_generator((self.commitments[0].len() as u32).trailing_zeros() as u64).unwrap();

        // Simulate the "fold" stage
        for root in self.commitments.iter() {
	    transcript.append_message(b"merkle_root", root.as_slice());
            let random_bytes = transcript.random();
            let random_element = F::from_digest(&random_bytes);
            random_elements.push(random_element);
            current_gen *= current_gen;
        }
	// Last fold step
        transcript.append_message(b"last_element", self.last_elem.as_ref());

        // Simulate the "query" stage
        for query in &self.queries {
            query.verify(
                &self.commitments,
                current_gen,
                self.last_elem,
                &random_elements,
            )?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ntt::NttField;
    use bincode; // Ensure bincode is imported

    #[test]
    fn fold_step_test() {
        let log_n = 5;
        let values: Vec<Field128> = (0..1 << log_n)
            .map(|i| Field128::from(i as i64 * 7 + 3))
            .collect();

        let gen = Field128::pow_2_generator(log_n + 1).unwrap();

        let mut transcript1 = Transcript::new();
        let mut transcript2 = Transcript::new();

        let mut prover1 = ProverData::init(values.clone(), gen, &mut transcript1);
        let mut prover2 = ProverData::init(values.clone(), gen, &mut transcript2);

        prover1.fold_step(gen, &mut transcript1);
        prover2.fold_step_opt(gen, &mut transcript2);

        assert_eq!(
            prover1.commitments[1].root(),
            prover2.commitments[1].root(),
            "Merkle roots differ after folding"
        );

        assert_eq!(
            prover1.commitments[1].data, prover2.commitments[1].data,
            "Commitment data differs after folding"
        );
    }

    // create a test for prove and verify!

    // create a test for creating the proof for a big RS code (1 million field
    // elements, or 16mb), serializes it and prints the size of the proof to
    // serialize please use `Serde` and `Bincode`. You'll have to add derive
    // instances for Serde in the proof datatype
    #[test]
    fn prove_and_verify_test() {
        let log_n = 5;
        let values: Vec<Field128> = (0..1 << log_n)
            .map(|i| Field128::from(i as i64 * 7 + 3))
            .collect();

        let mut transcript = Transcript::new();
        let proof = FriProof::prove(values.clone(), &mut transcript);

        assert!(proof.verify().is_ok(), "Proof verification failed");
    }

    #[test]
    fn big_rs_code_proof_test() {
        // Create a large RS code with 1 million elements
        let values: Vec<Field128> = (0..1 << 20).map(|i| Field128::from(i as i64)).collect();

        let mut transcript = Transcript::new();
        let proof = FriProof::prove(values.clone(), &mut transcript);

        // Serialize the proof using Serde and Bincode
        let serialized_proof = bincode::serialize(&proof).expect("Serialization failed");
        println!("Size of serialized proof: {} bytes", serialized_proof.len());

        // Deserialize the proof
        let _: FriProof<Field128> =
            bincode::deserialize(&serialized_proof).expect("Deserialization failed");
    }
}
