//! The crate's error type.

use std::fmt;

/// Why a modeling operation produced no result.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// The input violates the operation's preconditions (negative radius, unknown id,
    /// wrong kind of shape, ...).
    InvalidInput(String),
    /// The input was acceptable but the algorithm could not build a result.
    Failed(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            Error::Failed(msg) => write!(f, "operation failed: {msg}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;
