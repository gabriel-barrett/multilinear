use crate::field::Field;
use crate::field::Field128;
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
