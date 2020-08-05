//! Expanding Polytope Algorithm

pub use self::epa2d::{EPALeft2, EPA2};
pub use self::epa3d::{EPALeft3, EPAFn3, EPA3};

mod epa2d;
mod epa3d;

use cgmath::prelude::*;
use cgmath::{Point3, Vector3, BaseFloat};

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

    /// Create a new EPA instance
    fn new() -> Self;

    /// Create a new EPA instance with the given tolerance
    fn new_with_tolerance(
        tolerance: <Self::Point as EuclideanSpace>::Scalar,
        max_iterations: u32,
    ) -> Self;
}

/// A trait for a function that can be used to determine
/// which direction to resolve a collision in
pub trait EPAResolveDir3Fn<S: BaseFloat> {

    /// Gets the normal to report,
    /// given the initial normal calculated by EPA
    fn normal<SL, SR, TL, TR>(
        normal: &Vector3<S>,
        left: &SL,
        left_transform: &TL,
        right: &SR,
        right_transform: &TR,
    ) -> Vector3<S>
    where
        SL: Primitive<Point = Point3<S>>,
        SR: Primitive<Point = Point3<S>>,
        TL: Transform<Point3<S>>,
        TR: Transform<Point3<S>>;

    /// Gets the direction to resolve a collision in,
    /// given the normal to report.
    fn resolve_dir<SL, SR, TL, TR>(
        normal: &Vector3<S>,
        left: &SL,
        left_transform: &TL,
        right: &SR,
        right_transform: &TR,
    ) -> Vector3<S>
    where
        SL: Primitive<Point = Point3<S>>,
        SR: Primitive<Point = Point3<S>>,
        TL: Transform<Point3<S>>,
        TR: Transform<Point3<S>>;
}

/// A helper macro for defining a function for EPAFn3.
/// The syntax looks like this:
/// 
/// epa_resolve_dir3! {
///     pub EPALeftResolveDir<f64> = for<SL, TL, SR, TR>
///         normal |n, sl, tl, sr, tr| l.closest_valid_normal(n, lt),
///         resolve_dir |n, sl, tl, sr, tr| *n
/// }
/// 
#[macro_export]
macro_rules! epa_resolve_dir3 {
    {
        pub $name:ident < $S:ty > = for < $SL:ident, $TL:ident, $SR:ident, $TR:ident >
            normal |$n:ident, $l:ident, $lt:ident, $r:ident, $rt:ident| $f:expr,
            resolve_dir |$n2:ident, $l2:ident, $lt2:ident, $r2:ident, $rt2:ident| $f2:expr
    } => {
        epa_resolve_dir3! {
            (pub) $name<$S> = for<$SL, $TL, $SR, $TR>
                normal |$n, $l, $lt, $r, $rt| $f,
                resolve_dir |$n2, $l2, $lt2, $r2, $rt2| $f2
        };
    };

    {
        $name:ident < $S:ty > = for < $SL:ident, $TL:ident, $SR:ident, $TR:ident >
            normal |$n:ident, $l:ident, $lt:ident, $r:ident, $rt:ident| $f:expr,
            resolve_dir |$n2:ident, $l2:ident, $lt2:ident, $r2:ident, $rt2:ident| $f2:expr
    } => {
        epa_resolve_dir3! {
            () $name<$S> = for<$SL, $TL, $SR, $TR>
                normal |$n, $l, $lt, $r, $rt| $f,
                resolve_dir |$n2, $l2, $lt2, $r2, $rt2| $f2
        };
    };

    {
        ($($vis:tt)*) $name:ident < $S:ty > = for < $SL:ident, $TL:ident, $SR:ident, $TR:ident >
            normal |$n:ident, $l:ident, $lt:ident, $r:ident, $rt:ident| $f:expr,
            resolve_dir |$n2:ident, $l2:ident, $lt2:ident, $r2:ident, $rt2:ident| $f2:expr
    } => {
        $($vis)* struct $name;

        impl $crate::algorithm::minkowski::EPAResolveDir3Fn<$S> for $name {
            #[allow(unused_qualifications)]
            fn normal<$SL, $SR, $TL, $TR>(
                $n: &cgmath::Vector3<$S>,
                $l: &$SL,
                $lt: &$TL,
                $r: &$SR,
                $rt: &$TR,
            ) -> cgmath::Vector3<$S>
            where
                $SL: $crate::prelude::Primitive<Point = cgmath::Point3<$S>>,
                $SR: $crate::prelude::Primitive<Point = cgmath::Point3<$S>>,
                $TL: cgmath::prelude::Transform<cgmath::Point3<$S>>,
                $TR: cgmath::prelude::Transform<cgmath::Point3<$S>>
            {
                $f
            }

            #[allow(unused_qualifications)]
            fn resolve_dir<$SL, $SR, $TL, $TR>(
                $n2: &cgmath::Vector3<$S>,
                $l2: &$SL,
                $lt2: &$TL,
                $r2: &$SR,
                $rt2: &$TR,
            ) -> cgmath::Vector3<$S>
            where
                $SL: $crate::prelude::Primitive<Point = cgmath::Point3<$S>>,
                $SR: $crate::prelude::Primitive<Point = cgmath::Point3<$S>>,
                $TL: cgmath::prelude::Transform<cgmath::Point3<$S>>,
                $TR: cgmath::prelude::Transform<cgmath::Point3<$S>>
            {
                $f2
            }
        }
    };
}