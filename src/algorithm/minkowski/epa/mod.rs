//! Expanding Polytope Algorithm

pub use self::epa2d::{EPALeft2, EPA2};
pub use self::epa3d::{EPALeft3, EPA3};

mod epa2d;
mod epa3d;

use cgmath::prelude::*;

use super::SupportPoint;
use crate::prelude::*;
use crate::Contact;

pub const EPA_TOLERANCE: f32 = 0.00001;
pub const MAX_ITERATIONS: u32 = 100;

/// Expanding Polytope Algorithm base trait
pub trait EPA {
    /// Point type
    type Point: EuclideanSpace;

    /// Process the given simplex, and compute the contact point.
    ///
    /// The given simplex must be a complete simplex for the given space, and it must enclose the
    /// origin.
    fn process<SL, SR, TL, TR>(
        &self,
        simplex: &mut Vec<SupportPoint<Self::Point>>,
        left: &SL,
        left_transform: &TL,
        right: &SR,
        right_transform: &TR,
    ) -> Option<Contact<Self::Point>>
    where
        SL: Primitive<Point = Self::Point>,
        SR: Primitive<Point = Self::Point>,
        TL: Transform<Self::Point>,
        TR: Transform<Self::Point>;

    /// Process the given simplex, and compute the contact point.
    ///
    /// Takes a function that outputs the true normal given the input normal and shapes,
    /// and a function that outputs the resolve_dir given the normal and shapes.
    ///
    /// The given simplex must be a complete simplex for the given space, and it must enclose the
    /// origin.
    fn process_ex<SL, SR, TL, TR>(
        &self,
        simplex: &mut Vec<SupportPoint<Self::Point>>,
        left: &SL,
        left_transform: &TL,
        right: &SR,
        right_transform: &TR,
        _normal_fn: impl FnOnce(&<Self::Point as EuclideanSpace>::Diff, &SL, &TL, &SR, &TR) -> <Self::Point as EuclideanSpace>::Diff,
        _resolve_dir_fn: impl FnOnce(&<Self::Point as EuclideanSpace>::Diff, &SL, &TL, &SR, &TR) -> <Self::Point as EuclideanSpace>::Diff,
    ) -> Option<Contact<Self::Point>>
    where
        SL: Primitive<Point = Self::Point>,
        SR: Primitive<Point = Self::Point>,
        TL: Transform<Self::Point>,
        TR: Transform<Self::Point>
    {
        self.process(simplex, left, left_transform, right, right_transform)
    }

    /// Create a new EPA instance
    fn new() -> Self;

    /// Create a new EPA instance with the given tolerance
    fn new_with_tolerance(
        tolerance: <Self::Point as EuclideanSpace>::Scalar,
        max_iterations: u32,
    ) -> Self;
}