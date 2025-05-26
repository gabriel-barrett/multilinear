use serde::{Deserialize, Serialize};
use sha2::digest::{generic_array::GenericArray, OutputSizeUser};
use sha2::{Digest, Sha256}; // Import Serialize and Deserialize

pub type HashDigest = GenericArray<u8, <Sha256 as OutputSizeUser>::OutputSize>;

#[derive(Debug)]
pub struct Merkle<T> {
    pub layers: Vec<Vec<HashDigest>>,
    pub data: Vec<T>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Direction {
    Left,
    Right,
}

#[derive(Debug, Clone, Serialize, Deserialize)] // Add Serialize and Deserialize
pub struct MerkleInclusionPath<T> {
    pub value: T,
    pub path: Vec<(HashDigest, Direction)>,
}

impl<T> Merkle<T>
where
    T: AsRef<[u8]>,
{
    pub fn commit(data: Vec<T>) -> Self {
        assert!(
            data.len().is_power_of_two(),
            "Data length must be a power of two"
        );

        let first_layer: Vec<HashDigest> = data.iter().map(|item| hash_leaf(item)).collect();

        let mut layers: Vec<Vec<HashDigest>> = vec![first_layer];

        while layers.last().unwrap().len() > 1 {
            let current_layer = layers.last().unwrap();
            let next_layer: Vec<HashDigest> = current_layer
                .chunks(2)
                .map(|pair| hash_node(&pair[0], &pair[1]))
                .collect();
            layers.push(next_layer);
        }

        Merkle { layers, data }
    }

    pub fn root(&self) -> HashDigest {
        self.layers.last().unwrap()[0]
    }

    pub fn open(&self, index: usize) -> Option<MerkleInclusionPath<T>>
    where
        T: Clone,
    {
        if index >= self.data.len() {
            return None;
        }

        let value = self.data[index].clone();
        let mut path = Vec::new();
        let mut current_index = index;

        for layer in &self.layers {
            let (sibling_index, direction) = if current_index % 2 == 0 {
                (current_index + 1, Direction::Right)
            } else {
                (current_index - 1, Direction::Left)
            };

            if sibling_index < layer.len() {
                path.push((layer[sibling_index], direction));
            }

            current_index /= 2;
        }

        Some(MerkleInclusionPath { value, path })
    }
}

fn hash_leaf<T: AsRef<[u8]>>(item: &T) -> HashDigest {
    let mut hasher = Sha256::new();
    hasher.update(item.as_ref());
    hasher.finalize()
}

fn hash_node(left: &HashDigest, right: &HashDigest) -> HashDigest {
    let mut hasher = Sha256::new();
    hasher.update(left);
    hasher.update(right);
    hasher.finalize()
}

pub enum MerkleInclusionPathError {
    IncompatibleHash(HashDigest, HashDigest),
}

impl std::fmt::Debug for MerkleInclusionPathError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MerkleInclusionPathError::IncompatibleHash(hash1, hash2) => {
                write!(f, "{hash1:x?} != {hash2:x?}")
            }
        }
    }
}

impl<T> MerkleInclusionPath<T>
where
    T: AsRef<[u8]> + Clone,
{
    pub fn verify(&self, root: &HashDigest) -> Result<(), MerkleInclusionPathError> {
        let mut computed_hash = {
            let mut hasher = Sha256::new();
            hasher.update(self.value.as_ref());
            hasher.finalize()
        };

        for (sibling_hash, direction) in &self.path {
            computed_hash = match direction {
                Direction::Left => hash_node(sibling_hash, &computed_hash),
                Direction::Right => hash_node(&computed_hash, sibling_hash),
            };
        }

        if &computed_hash == root {
            Ok(())
        } else {
            Err(MerkleInclusionPathError::IncompatibleHash(
                *root,
                computed_hash,
            ))
        }
    }
}

#[cfg(test)]
#[test]
fn test_open_verify() {
    let data = vec![[0], [8], [4], [1], [5], [7], [6], [1]];
    let merkle_tree = Merkle::commit(data);

    println!("Merkle Root:\n {:x?}", merkle_tree.root());
    let proof = merkle_tree.open(5).unwrap();
    println!("Inclusion Path for index 5:\n {proof:x?}");
    proof.verify(&merkle_tree.root()).unwrap();
}
