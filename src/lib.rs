#![crate_type = "rlib"]
#![crate_type = "dylib"]
#![deny(
    missing_docs,
    trivial_casts,
    unstable_features,
    unused_import_braces,
    unused_qualifications
)]

//! Companion library to cgmath, dealing with collision detection centric data structures and
//! algorithms.
//!
//! This crate provides useful data structures and algorithms for doing collision detection.
//! It is organized into a few distinct parts: generic geometry (ray, line, plane, frustum etc),
//! bounding volumes (AABB, OBB, Sphere etc), collision primitives and algorithms used for
//! collision detection, distance computation etc.
//!

extern crate bit_set;
extern crate cgmath;
extern crate rand;

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
#[cfg_attr(test, macro_use)]
extern crate smallvec;
use wasm_bindgen::prelude::*;

// Re-exports

pub use bound::*;
pub use contact::*;
pub use frustum::*;
pub use line::*;
pub use plane::Plane;
pub use ray::*;
pub use traits::*;
pub use volume::*;

pub mod algorithm;
pub mod dbvt;
pub mod prelude;
pub mod primitive;

// Modules

mod bound;
mod contact;
mod frustum;
mod line;
mod plane;
mod ray;
mod traits;
mod volume;

/// Logs to the console.
#[macro_export]
macro_rules! console_log {
    ($($t:tt)*) => {
        $crate::log(&format_args!($($t)*).to_string())
    };
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}