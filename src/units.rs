#![warn(missing_docs)]
#![warn(clippy::pedantic)]

//! Crate for astronomy-specific unit support.  
//!
//! Many unit crates exist in rust and other programming languages, however, astronomy often
//! requires more esoteric units (solar mass, Mpc) and other specialized unit functionality on top
//! of general units. This crate aims to deliver a complete unit solution for astronomy.
//!
//! This package borrows from the [astropy.units](https://docs.astropy.org/en/stable/units/index.html) syntax so that users familiar with it should already
//! feel comfortable with how the units are implemented here. The syntax is very similar to how
//! units are actually written on paper.
//!
//! The crate comes with most units built in and building a quantity (a f64 value with a unit) can
//! be done through multiplication and division.  
//! ```
//! use astroxide::units::*;
//!
//! let distance = 5. * MEGAPARSEC;
//! assert_eq!(distance.value, 5.);
//! assert_eq!(distance.unit, MEGAPARSEC.as_unit());
//! ```
//!
//! Units of the same dimensionality can be added and subtracted and converted to one another
//! without worrying about the specific implementation.
//!
//! ```
//! use astroxide::units::*;
//! let speed_1 = 50. * KILOMETER/HOUR;
//! let speed_2 = 10. * METER/SECOND;
//! let sum = speed_1 + speed_2;
//! let sum_mph = sum.to(MILE/HOUR);
//! assert_eq!(sum.value, 86.);
//! assert!((sum_mph.value -53.4379).abs() < 0.001)
//! ```
//!
//! Units of differing dimensionality will panic if added, subtracted, or converted.
//!
//!
//! What makes this units package novel is the handling of cosmological units. That is to say,
//! taking into account the often dreaded [Little h](https://arxiv.org/pdf/1308.4150). astroxide
//! units allow us to factor out little h from existing quantities, assume values for little h, or
//! change cosmologies completely (In addition to all the standard conversions and dimension
//! checking that come with standard units).  
//!
//! For example if a paper assumes H<sub>0</sub> = 70 and quotes a distance measurement of 1 Mpc, then you can
//! factor out the little h dependency (assuming a standard h<sup>-1</sup> scaling for
//! distance) which on paper would be 0.7 h<sup>-1</sup> Mpc [see this article](https://www.astro.ljmu.ac.uk/~ikb/research/h-units.html).
//! This factored out value is a ``CosmoQuantity`` which is built of a ``CosmoValue`` (containing the little h dependency) and its exponent (-1 in this case).
//!
//! ```
//! use astroxide::units::*;
//!
//! let distance = 1. * MEGAPARSEC;
//! let factored_out_distance = distance.factor_out_h(0.7, -1);
//!
//! assert_eq!(factored_out_distance.cosmo_value.value, 0.7);
//! assert_eq!(factored_out_distance.cosmo_value.h_dependency, -1);
//! assert_eq!(factored_out_distance.unit, MEGAPARSEC.as_unit())
//! ```
//! Alternatively, if a paper instead has the dependency explicitly factored out, then we can do
//! the reverse and adopt a specific H<sub>0</sub> value. For example if we have a mass derived
//! from luminousity then the h scaling would go ~ h<sup>-2</sup>. We provide the little h struct
//! as syntactic sugar for including little h scalings.
//!
//! ```
//! use astroxide::units::*;
//!
//! let mass = 10.* h(-2) *SOLAR_MASS;
//! let mass_70 = mass.factor_in_h(0.7); // assuming H0 = 70 km/s/Mpc
//! assert!((mass_70.value- 20.4081632).abs() < 0.0001);
//! assert_eq!(mass_70.unit, SOLAR_MASS.as_unit());
//! ```
//!
//! Of course we can just change little h assumptions completely. If we have a volume (h<sup>-3</sup> scaling) assuming
//! H<sub>0</sub> = 70, and want to convert that to a cosmology assuming H<sub>0</sub> = 67, then
//! one would usually factor out the h dependency then assume a different h (which we can do above)
//! or simply use the ``Quantity::switch_cosmologies`` method on any quantity.
//!
//! ```
//! use astroxide::units::*;
//!
//! let volume_70 = 2000.* MEGAPARSEC * MEGAPARSEC * MEGAPARSEC;
//! let volume_67 = volume_70.switch_cosmologies(0.7, 0.67, -3);
//! assert!((volume_67.value - 2280.86566499).abs() < 0.0001);
//! ```
//!
//!
//!
use colored::Colorize;
use fmtastic::Superscript;
use paste::paste;
use std::cmp::Ordering;
use std::{
    collections::BTreeMap,
    fmt::{self, Display},
    ops::{Add, Div, Mul, Neg, Sub},
};

/// The dimensional "finger print" which defines the dimensionality of a unit.
///
/// This struct fully characterizes the dimensionality of unit which is essential for making sure
/// that units of the same dimensionality can be added or subtracted or converted to one another.
///
/// For example a speed of 5 km/hour and a speed of 10 miles/hour can be added togther and coverted
/// to one another because they both have a dimensionality of L<sup>1</sup> T<sup>-1</sup> (Length
/// over time). The `Dimension` Struct captures this.
///
/// It is worth noting that "angular distance" and "solid angle" are not, strictly speaking,
/// dimensions. However, for astronomical purposes which deal frequently with spherical
/// trigonometry, it is useful to treat them as such.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Dimension {
    /// scaling of "Length"" dimension.
    pub length: i32,
    /// scaling of "Mass" dimension.
    pub mass: i32,
    /// scaling of "Time" dimension.
    pub time: i32,
    /// scaling of "Temperature" dimension.
    pub temperature: i32,
    /// scaling of "Current" dimension.
    pub current: i32,
    /// scaling of "angular distance" dimension.
    pub angular_distance: i32,
    /// scaling of "solid angle" dimension.
    pub solid_angle: i32,
    /// scaling of "Luminous Intensity" dimension.
    pub luminous_intensity: i32,
    /// scaling of "amount of substance" dimension.
    pub amount_of_substance: i32,
}

impl Dimension {
    /// Helper method for creating a Dimensionless Dimension.
    pub const ZERO: Dimension = Dimension {
        length: 0,
        mass: 0,
        time: 0,
        temperature: 0,
        current: 0,
        angular_distance: 0,
        solid_angle: 0,
        luminous_intensity: 0,
        amount_of_substance: 0,
    };
    /// New method which creates an empty dimension.
    #[must_use]
    pub fn new() -> Self {
        Self::ZERO
    }
}

impl Default for Dimension {
    fn default() -> Self {
        Dimension::new()
    }
}

impl Mul<i32> for Dimension {
    type Output = Self;
    fn mul(self, rhs: i32) -> Self::Output {
        Dimension {
            length: self.length * rhs,
            mass: self.mass * rhs,
            time: self.time * rhs,
            temperature: self.temperature * rhs,
            current: self.current * rhs,
            angular_distance: self.angular_distance * rhs,
            solid_angle: self.solid_angle * rhs,
            luminous_intensity: self.luminous_intensity * rhs,
            amount_of_substance: self.amount_of_substance * rhs,
        }
    }
}
impl Add for Dimension {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Dimension {
            length: self.length + rhs.length,
            mass: self.mass + rhs.mass,
            time: self.time + rhs.time,
            temperature: self.temperature + rhs.temperature,
            current: self.current + rhs.current,
            angular_distance: self.angular_distance + rhs.angular_distance,
            solid_angle: self.solid_angle + rhs.solid_angle,
            luminous_intensity: self.luminous_intensity + rhs.luminous_intensity,
            amount_of_substance: self.amount_of_substance + rhs.amount_of_substance,
        }
    }
}

impl Sub for Dimension {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Dimension {
            length: self.length - rhs.length,
            mass: self.mass - rhs.mass,
            time: self.time - rhs.time,
            temperature: self.temperature - rhs.temperature,
            current: self.current - rhs.current,
            angular_distance: self.angular_distance - rhs.angular_distance,
            solid_angle: self.solid_angle - rhs.solid_angle,
            luminous_intensity: self.luminous_intensity - rhs.luminous_intensity,
            amount_of_substance: self.amount_of_substance - rhs.amount_of_substance,
        }
    }
}

/// Building block unit which can make more complex units.
///
/// Base units represent the indivisible units which make up other units and include the symbol for
/// the units, their dimensionality and their conversion factors to some arbitary unit within that
/// dimension.
///
/// For example. A unit like km/s is made up of two base units: kilometer and second which have the
/// symbols "km" and "s" with dimensions L and T respectively. Other units which can be made up of
/// other units can also be considered base units. A base unit for energy, for example, would be
/// Joule which has the symbol "J" and the dimensions of M L<sup>2</sup> T<sup>-2</sup>. But an
/// equivalent unit is kg m<sup>2</sup> s<sup>-2</sup>.
///
/// So base units are specifically modeling the unit representations which can be combined with other units to make new units and NOT units with a single dimensionality.
#[derive(Debug, Clone, Copy)]
pub struct BaseUnit {
    /// Dimensional finger print of unit.
    pub base_dimension: Dimension,
    /// The symbol representation of the unit e.g., "km", "s", "Mpc".
    pub symbol: &'static str,
    /// The conversion factor to some arbitary common value for that dimensionality. Usually the SI standard.
    pub conversion_factor: f64,
}

impl Ord for BaseUnit {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.base_dimension, &self.symbol).cmp(&(other.base_dimension, &other.symbol))
    }
}

impl PartialOrd for BaseUnit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for BaseUnit {
    fn eq(&self, other: &Self) -> bool {
        (self.base_dimension == other.base_dimension) && (self.symbol == other.symbol)
    }
}

impl Eq for BaseUnit {}

/// Unit represntation.
///
/// Units are represented as a simple ``std::collections::BTreeMap`` of Base Units and their exponents. This allows for
/// any arbitary unit to be manually defined from the base units that exist. For example, square
/// meters are simply (meters: 2) whereas mpc would be defined as (mile: 1, hour: -1).
///
/// This allows us to compare if units are equal as well as quickly derived the dimensionality of a
/// unit, taking into account the exponents of the base units. So we can say that square feet and
/// square meters are not the same unit (foot: 2) != (meter: 2), do have the same dimensionality.
#[derive(Debug, Clone)]
pub struct Unit {
    /// Mapping of Base units and their corresponding exponents.
    /// e.g. km/s would be `BTreeMap::from([(KILOMETER, 1), (SECOND, -1)])` which would be
    /// equivalent to `KILOMETER/SECOND`.
    /// ```
    /// use astroxide::units::*;
    /// use std::collections::BTreeMap;
    /// let unit_a = Unit {base_units: BTreeMap::from([(KILOMETER, 1), (SECOND, -1)])};
    /// let unit_b = KILOMETER/SECOND;
    /// assert_eq!(unit_a, unit_b);
    /// ```
    pub base_units: BTreeMap<BaseUnit, i32>,
}

/// Main data structure for representing values and their associated units.
///
/// Quantity combines f64 values and their units into one struct. This is the main way that a user
/// would interact with data. The 'quantity' 5 km is distinct from the specific value '5' and the abstract unit 'km'.
#[derive(Debug, Clone, PartialEq)]
pub struct Quantity {
    /// The unit of the quantity represented as [Unit] struct. ('km/s' in the quantity 5 km/s).
    pub unit: Unit,
    /// f64 value of the quantity. ('5' in the quantity 5 km/s)
    pub value: f64,
}

/// Describes things that can act like units.
///
/// Things that are directly convertable to units can use this trait. This is mainly used for
/// ``BaseUnits`` to be treated as regular units. This also includes the calculation of the converson
/// factor for the thing implementing the trait.
pub trait UnitLike {
    /// Conversion to a Unit which makes the thing "unitlike"
    fn as_unit(&self) -> Unit;
    /// Default implementation which calculates the conversion factor. For example if we have
    /// km<sup>2</sup> as a Unit (which trivially implements Unitlike) then the conversion factor for
    /// the unit would be the conversion factor of the base unit(km) to the power of the exponent
    /// of the unit (2).
    fn calculate_conversion_factor(&self) -> f64 {
        self.as_unit()
            .base_units
            .iter()
            .map(|(base_unit, &exponent)| base_unit.conversion_factor.powi(exponent))
            .product()
    }
}

impl Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, (base_unit, &exponent)) in self.base_units.iter().enumerate() {
            if i > 0 {
                write!(f, " ")?;
            }
            match exponent {
                1 => write!(f, "{}", base_unit.symbol)?,
                _ => write!(f, "{}{}", base_unit.symbol, Superscript(exponent))?,
            }
        }
        Ok(())
    }
}

impl Unit {
    /// The dimensionality of the given unit.
    ///
    /// Every unit can be defined by its dimensonal "fingerprint". This method builds this
    /// fingerprint from the base unit dimensions which build up the unit. This is critical for
    /// ensuring that units of the same dimension can be added, subtracted, and converted.
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let x = METER/SECOND;
    /// let dimensions = x.dimensions();
    /// assert_eq!(dimensions, dim!(length: 1, time: -1));
    /// ```
    #[must_use]
    pub fn dimensions(&self) -> Dimension {
        let mut full_dimension = Dimension::new();
        for (base_unit, &exponent) in &self.base_units {
            full_dimension = full_dimension + (base_unit.base_dimension * exponent);
        }
        full_dimension
    }
    /// Inverse of the Unit.
    ///
    /// This is equivalent to return the unit to the power of -1.
    /// For example, if we have the unit "km/s" then the inverse would return "s/km".
    /// ```
    /// use astroxide::units::*;
    ///
    /// let x = METER/SECOND;
    /// let x_inverse = x.invert();
    /// assert_eq!(x_inverse, SECOND/METER);
    /// ```
    /// This is helpful for implementing divisions.
    ///
    #[must_use]
    pub fn invert(&self) -> Self {
        let mut inverted = BTreeMap::new();
        for (base_unit, exponent) in &self.base_units {
            inverted.insert(*base_unit, -exponent);
        }
        Unit {
            base_units: inverted,
        }
    }
}

impl Quantity {
    /// Creates a new Quantity object.
    ///
    /// This is the explicit way to create a Quantity object which can be useful but is fully
    /// equivalent to the syntactic sugar of multiplying the unit by a f64.
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let unit_1 = 5. * METER;
    /// let unit_2 = Quantity::new(5., METER);
    /// assert_eq!(unit_1, unit_2);
    /// ```
    pub fn new(value: f64, unit: impl UnitLike) -> Self {
        Self {
            value,
            unit: unit.as_unit(),
        }
    }
    /// Conversion from one unit to an equivalent unit.
    ///
    /// Core functionality of [Quantity] is to allow conversion from one unit to an equivalent
    /// one. This is done using the `to` method.
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let speed = 5.* METER/SECOND;
    /// let speed_mph = speed.to(MILE/HOUR);
    /// let answer = 11.1847*(MILE/HOUR);
    /// assert!(speed_mph.approx_eq(&answer, 3));
    /// ```
    /// # Panics
    /// `.to` will result in a panic if the units do not have the dame dimensionality.
    /// ```should_panic(expected = Cannot convert)
    /// use astroxide::units::*;
    /// let speed = 5. *MEGAPARSEC/HOUR;
    /// let area = GIGAPARSEC*GIGAPARSEC;
    /// let conversion = speed.to(area); // non-physical conversion.
    /// ```
    #[must_use]
    pub fn to(&self, target_unit: impl UnitLike) -> Self {
        assert!(
            self.unit.dimensions() == target_unit.as_unit().dimensions(),
            "Cannot convert {} to {} since they have different dimensions.",
            self.unit,
            target_unit.as_unit()
        );

        Quantity {
            unit: target_unit.as_unit(),
            value: self.value
                * (self.unit.calculate_conversion_factor()
                    / target_unit.calculate_conversion_factor()),
        }
    }
    /// Approximately equal method for a given decimal precision.
    ///
    /// Helper method for determining if two quantities are approximately equal. Since units are
    /// only using f64 floating point errors can occur. In addition there are some use cases where
    /// absolute precision is not possible (comparison to other quantities in other works with
    /// different significant figures for example.), Here `decimal_place` refers to up to what
    /// decimal should be included in the comparison.
    ///
    /// Units are included in the comparison, so even if the values are within the tolerance,
    /// non-equal units will result in a "false".
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let x = 5.* METER;
    /// let y = 5.001 * METER;
    /// assert!(x.approx_eq(&y, 2));
    ///
    /// let unit = 5. * METER;
    /// let different_unit = 5. * KILOMETER;
    /// assert!(!unit.approx_eq(&different_unit, 3));
    /// ```
    #[must_use]
    pub fn approx_eq(&self, other: &Self, decimal_place: i32) -> bool {
        ((self.value / other.value) - 1.).abs() < 0.1_f64.powi(decimal_place)
            && self.unit == other.unit
    }

    /// Removes the assumed *h* dependency.
    ///
    /// *Opposite of [`CosmoQuantity::factor_in_h`].*
    ///
    /// Some quantities have been derived because of explicit assumptions of the cosmology that was
    /// used in the measurement of that quantity. For example, a Mass derived from Luminosity might
    /// well be determined to be 1.5e15 M<sub>☉</sub>. However, this quantity, without an explicit
    /// *h* scaling, usually means: 1. a *h*<sup>-2</sup> dependency and 2. a specific choice of *h*.
    ///
    /// In this example, if we knew that the assumed value of *h* was 0.7, then we could "factor
    /// out" *h* and write the mass as the *h*-independent  7.35e14 *h*<sup>-2</sup> M<sub>☉</sub>.
    /// This now shows the *h* scaling explicitly and allows another assumption of *h* to be used. (see [(Croton et. al., 2013)](https://arxiv.org/abs/1308.4150) and in particular "Case 4" therein.)
    ///
    /// In order to represent this quantity independent of any specific assumption of *h* we can
    /// use the `factor_out_h` function to explicit factor out the h assumption [(see examples here)](voluastro.ljmu.ac.uk/~ikb/research/h-units.html)
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let mass = 1.5e15*SOLAR_MASS;
    /// let mass_h_independent = mass.factor_out_h(0.7, -2);
    /// let answer = 7.35e14 * h(-2) * SOLAR_MASS;
    /// assert!(mass_h_independent.approx_eq(&answer,2));
    ///
    /// // examples from https://www.astro.ljmu.ac.uk/~ikb/research/h-units.html
    /// let value = 1. * MEGAPARSEC;
    /// let value_h = value.factor_out_h(0.7, -1);
    /// let answer = 0.7 * h(-1) * MEGAPARSEC;
    /// assert!(value_h.approx_eq(&answer, 9));
    ///
    /// ```
    #[must_use]
    pub fn factor_out_h(&self, h_value: f64, h_dependency: i32) -> CosmoQuantity {
        CosmoQuantity::new(
            self.value / (h_value.powi(h_dependency)),
            h_dependency,
            self.unit.clone(),
        )
    }
    /// Converts the quantity to a different assumption of *h*.
    ///
    /// It is often necessary to convert quanities in one assumped cosmology into another. For
    /// example, a distance of 1 Mpc assuming *h* = 0.7 and a scaling of *h<sup>-1</sup>* can be
    /// converted into a cosmology where *h* = 0.67 by first factoring out the *h* dependency and
    /// then factoring in a different assumption for *h*. This is often essential for comparing
    /// identical physical properties that might vary in *h* scaling (Mass from luminosity vs mass
    /// from velocity dispersion, for example).
    ///
    /// Switch cosmologies is a helper function which allows changing `from_h` `to_h` for a given h
    /// dependency.
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let distance_70 = 1.*MEGAPARSEC;
    /// let distance_60 = distance_70.switch_cosmologies(0.7, 0.6, -1);
    /// let answer = 1.1667 * MEGAPARSEC;
    /// assert!(distance_60.approx_eq(&answer, 2))
    /// ```
    #[must_use]
    pub fn switch_cosmologies(&self, from_h: f64, to_h: f64, h_dependency: i32) -> Self {
        let with_h = self.factor_out_h(from_h, h_dependency);
        with_h.factor_in_h(to_h)
    }
}

impl Display for Quantity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.unit)
    }
}

impl UnitLike for Unit {
    fn as_unit(&self) -> Unit {
        self.to_owned()
    }
}

impl UnitLike for BaseUnit {
    fn as_unit(&self) -> Unit {
        Unit {
            base_units: BTreeMap::from([(*self, 1)]),
        }
    }
}

impl PartialEq for Unit {
    fn eq(&self, other: &Self) -> bool {
        self.base_units == other.base_units
    }
}

impl Div<BaseUnit> for BaseUnit {
    type Output = Unit;
    fn div(self, rhs: BaseUnit) -> Self::Output {
        self.as_unit() / rhs.as_unit()
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div<Unit> for Unit {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.invert()
    }
}

impl Div<Quantity> for Quantity {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let unit = self.unit / rhs.unit;
        let value = self.value / rhs.value;

        Quantity { unit, value }
    }
}

impl Div<Unit> for BaseUnit {
    type Output = Unit;
    fn div(self, rhs: Unit) -> Self::Output {
        self.as_unit() / rhs
    }
}
impl Div<BaseUnit> for Unit {
    type Output = Unit;
    fn div(self, rhs: BaseUnit) -> Self::Output {
        self / rhs.as_unit()
    }
}

impl Div<f64> for BaseUnit {
    type Output = Quantity;
    fn div(self, rhs: f64) -> Self::Output {
        Quantity {
            value: 1. / rhs,
            unit: self.as_unit(),
        }
    }
}

impl Div<BaseUnit> for f64 {
    type Output = Quantity;
    fn div(self, rhs: BaseUnit) -> Self::Output {
        Quantity {
            value: self,
            unit: rhs.as_unit().invert(),
        }
    }
}

impl Div<f64> for Quantity {
    type Output = Quantity;
    fn div(self, rhs: f64) -> Self::Output {
        Quantity {
            value: self.value / rhs,
            unit: self.unit,
        }
    }
}
impl Div<Quantity> for f64 {
    type Output = Quantity;
    fn div(self, rhs: Quantity) -> Self::Output {
        Quantity {
            value: self / rhs.value,
            unit: rhs.unit.invert(),
        }
    }
}

impl Div<BaseUnit> for Quantity {
    type Output = Quantity;
    fn div(self, rhs: BaseUnit) -> Self::Output {
        let new_unit = self.unit / rhs.as_unit();
        Quantity {
            value: self.value,
            unit: new_unit,
        }
    }
}
impl Div<Quantity> for BaseUnit {
    type Output = Quantity;
    fn div(self, rhs: Quantity) -> Self::Output {
        let new_unit = self.as_unit() / rhs.unit;
        Quantity {
            value: 1. / rhs.value,
            unit: new_unit,
        }
    }
}

impl Mul for BaseUnit {
    type Output = Unit;
    fn mul(self, rhs: Self) -> Self::Output {
        self.as_unit() * rhs.as_unit()
    }
}

impl Mul for Unit {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut new_units = rhs.base_units;

        for (base_unit, exponent) in self.base_units {
            new_units
                .entry(base_unit)
                .and_modify(|curr| *curr += exponent)
                .or_insert(exponent);
        }

        new_units.retain(|_, exponent| *exponent != 0);
        Unit {
            base_units: new_units,
        }
    }
}

impl Mul for Quantity {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let unit = self.unit * rhs.unit;
        let value = self.value * rhs.value;
        Quantity { unit, value }
    }
}

impl Mul<f64> for BaseUnit {
    type Output = Quantity;
    fn mul(self, rhs: f64) -> Self::Output {
        let unit: Unit = Unit {
            base_units: BTreeMap::from([(self, 1)]),
        };
        Quantity { unit, value: rhs }
    }
}

impl Mul<BaseUnit> for f64 {
    type Output = Quantity;
    fn mul(self, rhs: BaseUnit) -> Self::Output {
        rhs * self
    }
}

impl Mul<BaseUnit> for Unit {
    type Output = Self;
    fn mul(self, rhs: BaseUnit) -> Self::Output {
        self * rhs.as_unit()
    }
}

impl Mul<Unit> for BaseUnit {
    type Output = Unit;
    fn mul(self, rhs: Unit) -> Self::Output {
        self.as_unit() * rhs
    }
}

impl Mul<f64> for Unit {
    type Output = Quantity;
    fn mul(self, rhs: f64) -> Self::Output {
        Quantity {
            value: rhs,
            unit: self,
        }
    }
}
impl Mul<Unit> for f64 {
    type Output = Quantity;
    fn mul(self, rhs: Unit) -> Self::Output {
        Quantity {
            value: self,
            unit: rhs,
        }
    }
}

impl Mul<BaseUnit> for Quantity {
    type Output = Quantity;
    fn mul(self, rhs: BaseUnit) -> Self::Output {
        let new_unit = self.unit * rhs.as_unit();
        Quantity {
            value: self.value,
            unit: new_unit,
        }
    }
}
impl Mul<Quantity> for BaseUnit {
    type Output = Quantity;
    fn mul(self, rhs: Quantity) -> Self::Output {
        rhs * self
    }
}
impl Mul<f64> for Quantity {
    type Output = Quantity;
    fn mul(self, rhs: f64) -> Self::Output {
        Quantity {
            value: rhs * self.value,
            unit: self.unit,
        }
    }
}
impl Mul<Quantity> for f64 {
    type Output = Quantity;
    fn mul(self, rhs: Quantity) -> Self::Output {
        rhs * self
    }
}

impl Add for Quantity {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        assert!(
            self.unit.dimensions() == rhs.unit.dimensions(),
            "Cannot add or subtract {} and {} since they have different dimensions.",
            self.unit,
            rhs.unit,
        );
        let new_value = (self.value * self.unit.calculate_conversion_factor()
            + rhs.value * rhs.unit.calculate_conversion_factor())
            / self.unit.calculate_conversion_factor();
        Quantity {
            unit: self.unit,
            value: new_value,
        }
    }
}

impl Neg for Quantity {
    type Output = Quantity;
    fn neg(self) -> Self::Output {
        Quantity {
            value: -self.value,
            unit: self.unit,
        }
    }
}

impl Sub for Quantity {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

/// Cosmological value that contains an explicit *h* scaling.
#[derive(Debug, Clone, PartialEq)]
pub struct CosmoValue {
    /// Floating value which is to be scaled by *h*.
    pub value: f64,
    /// Explicit h scaling power. For example, 1 h<sup>-1</sup> will have a `h_dependency` of -1.
    pub h_dependency: i32,
}

impl CosmoValue {
    /// Approximately equal method for a given decimal precision.
    ///
    /// Helper method for determining if two `CosmoValue` are approximately equal. Since units are
    /// only using f64 floating point errors can occur. In addition there are some use cases where
    /// absolute precision is not possible (comparison to other quantities in other works with
    /// different significant figures for example.), Here `decimal_place` refers to up to what
    /// decimal should be included in the comparison.
    ///
    /// *h* dependency is included in the comparison, so even if the values are within the tolerance,
    /// non-equivalent h dependency will result in a "false".
    ///
    /// ```
    /// use astroxide::units::*;
    /// let x = CosmoValue::new(5., -1);
    /// let y = CosmoValue::new(5., -2);
    /// let z = CosmoValue::new(5.0001, -1);
    /// assert!(!x.approx_eq(&y, 9));
    /// assert!(x.approx_eq(&z, 3));
    /// ```
    #[must_use]
    pub fn approx_eq(&self, other: &Self, decimal_place: i32) -> bool {
        ((self.value / other.value) - 1.).abs() < 0.1_f64.powi(decimal_place)
            && self.h_dependency == other.h_dependency
    }
}

impl<T: UnitLike> Mul<T> for CosmoValue {
    type Output = CosmoQuantity;
    fn mul(self, rhs: T) -> Self::Output {
        CosmoQuantity {
            cosmo_value: self,
            unit: rhs.as_unit(),
        }
    }
}

impl CosmoValue {
    /// Constructs a new `CosmoValue` from the given value and `h_dependency`.
    ///
    /// # Examples
    /// Creating the value 0.7 *h*<sup>-1</sup>. Can be done by:
    ///
    /// ```
    /// use astroxide::units::*;
    /// let x = CosmoValue::new(0.7, -1);
    /// ```
    #[must_use]
    pub fn new(value: f64, h_dependency: i32) -> Self {
        Self {
            value,
            h_dependency,
        }
    }
}

impl Display for CosmoValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {}{}",
            self.value,
            "h".italic(),
            Superscript(self.h_dependency)
        )
    }
}

/// Struct representation of little h (*h*).
///
/// This is pure syntactic sugar for building [`CosmoValue`] objects and calls [`CosmoValue::new`] in
/// the background.
///
/// # Examples
/// Creating 0.7 *h*<sup>-1</sup>:
/// ```
/// use astroxide::units::*;
///
/// let cosmo_val_sugar = 0.7 * h(-1);
/// let cosmo_val_normal =  CosmoValue::new(0.7, -1);
/// assert_eq!(cosmo_val_sugar, cosmo_val_sugar);
/// ```
/// This is especially nice when buidling [`CosmoQuantity`] objects which also contain a unit.
/// The quantitiy 4 *h*<sup>-1</sup> Mpc can be represented in a nearly identical way.
/// ```
/// use astroxide::units::*;
///
/// let cosmo_quantitiy = 4. *h(-1) * MEGAPARSEC;
/// assert_eq!(cosmo_quantitiy.cosmo_value, 4. * h(-1));
/// assert_eq!(cosmo_quantitiy.unit, MEGAPARSEC.as_unit());
/// ```
///
#[allow(non_camel_case_types)]
pub struct h(pub i32);

impl Mul<f64> for h {
    type Output = CosmoValue;
    fn mul(self, rhs: f64) -> Self::Output {
        CosmoValue::new(rhs, self.0)
    }
}

impl Mul<h> for f64 {
    type Output = CosmoValue;
    fn mul(self, rhs: h) -> Self::Output {
        CosmoValue::new(self, rhs.0)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Mul for h {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        h(self.0 + rhs.0)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Mul<CosmoValue> for h {
    type Output = CosmoValue;
    fn mul(self, rhs: CosmoValue) -> Self::Output {
        CosmoValue::new(rhs.value, rhs.h_dependency + self.0)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Mul<h> for CosmoValue {
    type Output = CosmoValue;
    fn mul(self, rhs: h) -> Self::Output {
        CosmoValue::new(self.value, self.h_dependency + rhs.0)
    }
}

impl Div<f64> for h {
    type Output = CosmoValue;
    fn div(self, rhs: f64) -> Self::Output {
        CosmoValue::new(1. / rhs, self.0)
    }
}

impl Div<h> for f64 {
    type Output = CosmoValue;
    fn div(self, rhs: h) -> Self::Output {
        CosmoValue::new(self, -rhs.0)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for h {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        h(self.0 - rhs.0)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div<h> for CosmoValue {
    type Output = CosmoValue;
    fn div(self, rhs: h) -> Self::Output {
        CosmoValue::new(self.value, self.h_dependency - rhs.0)
    }
}

impl Div<CosmoValue> for h {
    type Output = CosmoValue;
    fn div(self, rhs: CosmoValue) -> Self::Output {
        CosmoValue::new(1. / rhs.value, self.0 - rhs.h_dependency)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<T: UnitLike> Div<T> for CosmoValue {
    type Output = CosmoQuantity;
    fn div(self, rhs: T) -> Self::Output {
        self * rhs.as_unit().invert()
    }
}

/// Equivalent to [`Quantity`] except with a [`CosmoValue`] instead of a primitive f64.
///
/// This type represents quantities that have a value (f64), a h scaling (i32) and a unit ([`Unit`]).
/// For example, 1 *h*<sup>-1</sup> Mpc would be considered a [`CosmoQuantity`] where as 1 Mpc is
/// just a [`Quantity`].
#[derive(Debug, Clone, PartialEq)]
pub struct CosmoQuantity {
    /// The f64 value and its *h* scaling defined by a [`CosmoValue`].
    pub cosmo_value: CosmoValue,
    /// The unit of the `CosmoQuantity`.
    pub unit: Unit,
}

impl CosmoQuantity {
    /// Constructs a new `CosmoQuantity` from the given value, `h_dependency`, and unit.
    ///
    /// # Examples
    /// Creating the value 0.7 *h*<sup>-1</sup> Mpc. Can be done by:
    ///
    /// ```
    /// use astroxide::units::*;
    /// let x = CosmoQuantity::new(0.7, -1, MEGAPARSEC);
    /// ```
    pub fn new(value: f64, h_dependency: i32, unit: impl UnitLike) -> Self {
        Self {
            cosmo_value: CosmoValue {
                value,
                h_dependency,
            },
            unit: unit.as_unit(),
        }
    }
    /// Approximately equal method for a given decimal precision.
    ///
    /// Helper method for determining if two comso quantities are approximately equal. Since units are
    /// only using f64 floating point errors can occur. In addition there are some use cases where
    /// absolute precision is not possible (comparison to other quantities in other works with
    /// different significant figures for example.), Here `decimal_place` refers to up to what
    /// decimal should be included in the comparison.
    ///
    /// Units and h dependencies are included in the comparison, so even if the values are within the tolerance,
    /// non-equal units and dependencies will result in a "false".
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let x = CosmoQuantity::new(0.7, -1, MEGAPARSEC);
    /// let y = CosmoQuantity::new(0.7, -2, MEGAPARSEC);
    /// let z = CosmoQuantity::new(0.7, -1, GIGAPARSEC);
    ///
    /// let a = CosmoQuantity::new(0.700001, -1, MEGAPARSEC);
    ///
    /// assert!(!x.approx_eq(&y, 9)); // different h scaling.
    /// assert!(!x.approx_eq(&z, 9)); // different units.
    ///
    /// assert!(x.approx_eq(&a, 4));
    /// ```
    #[must_use]
    pub fn approx_eq(&self, other: &Self, decimal_place: i32) -> bool {
        self.cosmo_value
            .approx_eq(&other.cosmo_value, decimal_place)
            && self.unit == other.unit
    }
    /// Makes an explicit asumption for *h*.
    ///
    /// *Opposite of [`Quantity::factor_out_h`].*
    ///
    /// This method is essential a way to convert from [`CosmoQuantity`] to a regular [`Quantity`] by
    /// making an assumption for *h*. For example, 1 *h*<sup>-1</sup> Mpc has no asumption of *h*,
    /// however, we could adopt a value of *h*=0.7, which would result in a quantity of 1 * (0.7)<sup>-1</sup> Mpc = 1.4285 Mpc.
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let value_h = 1. * h(-1) * MEGAPARSEC;
    /// let value = value_h.factor_in_h(0.7);
    /// let answer = 1.4285 * MEGAPARSEC;
    /// assert!(value.approx_eq(&answer, 3));
    ///
    /// ```
    #[must_use]
    pub fn factor_in_h(&self, h_value: f64) -> Quantity {
        Quantity {
            unit: self.unit.clone(),
            value: self.cosmo_value.value * (h_value).powi(self.cosmo_value.h_dependency),
        }
    }

    // Ignores that h scaling and just returns the value and the unit.
    // Useful for doing the unit conversions but is not physical and so kept private to avoid
    // confusion.
    fn ignore_h(&self) -> Quantity {
        Quantity {
            unit: self.unit.clone(),
            value: self.cosmo_value.value,
        }
    }

    /// Inverse of the `CosmoQuantity`
    ///
    /// This is equivalent to return the `CosmoQuantity` to the power of -1.
    /// For example, if we have the cosmo quantity 2 *h*<sup>-1</sup> Mpc then
    /// the inverse would be 0.5 *h* Mpc<sup>-1</sup>.
    /// ```
    /// use astroxide::units::*;
    ///
    /// let x = 2. * h(-1) * MEGAPARSEC;
    /// let x_inverse = x.invert();
    /// let answer = 0.5 * h(1) /MEGAPARSEC ;
    /// assert_eq!(x_inverse, answer);
    /// ```
    /// This is helpful for implementing divisions.
    ///
    #[must_use]
    pub fn invert(&self) -> CosmoQuantity {
        CosmoQuantity {
            cosmo_value: CosmoValue {
                value: 1. / self.cosmo_value.value,
                h_dependency: -self.cosmo_value.h_dependency,
            },
            unit: self.unit.clone().invert(),
        }
    }
    /// Convert to equivalent unit (and h dependency).
    ///
    /// This is the [`CosmoQuantity`] equivalent of [`Quantity::to`].
    ///
    /// Just like a [`Quantity`], [`CosmoQuantity`] can be converted to an equivalent unit. For
    /// example, 1000 *h*<sup>-1</sup> Mpc is equivalent to 1 *h*<sup>-1</sup> Gpc. The *h* scaling
    /// carries and only the units and value are effected.
    ///
    /// ```
    /// use astroxide::units::*;
    ///
    /// let distance_mpc = 1000. * h(-1) *MEGAPARSEC;
    /// let distance_gpc = distance_mpc.to(GIGAPARSEC);
    /// let answer = 1. * h(-1) * GIGAPARSEC;
    ///
    /// assert!(distance_gpc.approx_eq(&answer, 9));
    /// ```
    #[must_use]
    pub fn to(&self, unit: impl UnitLike) -> CosmoQuantity {
        let converted = self.ignore_h().to(unit);
        CosmoQuantity::new(
            converted.value,
            self.cosmo_value.h_dependency,
            converted.unit,
        )
    }
}

impl Display for CosmoQuantity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.cosmo_value, self.unit)
    }
}

impl Neg for CosmoQuantity {
    type Output = Self;
    fn neg(self) -> Self::Output {
        CosmoQuantity {
            cosmo_value: CosmoValue {
                value: -self.cosmo_value.value,
                h_dependency: self.cosmo_value.h_dependency,
            },
            unit: self.unit,
        }
    }
}

impl Add for CosmoQuantity {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        assert!(
            self.cosmo_value.h_dependency == rhs.cosmo_value.h_dependency,
            "Cannot add or subtract {self} {rhs} because of differing h dependencies."
        );
        let quant = self.ignore_h() + rhs.ignore_h();
        CosmoQuantity {
            cosmo_value: CosmoValue {
                value: quant.value,
                h_dependency: self.cosmo_value.h_dependency,
            },
            unit: quant.unit,
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Sub for CosmoQuantity {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for CosmoQuantity {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        CosmoQuantity {
            cosmo_value: CosmoValue {
                value: self.cosmo_value.value * rhs.cosmo_value.value,
                h_dependency: self.cosmo_value.h_dependency + rhs.cosmo_value.h_dependency,
            },
            unit: self.unit * rhs.unit,
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for CosmoQuantity {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.invert()
    }
}

/// Macro for creating [Dimension] object.
///
/// Equivalent to creating a Dimension struct but since all fields need a value (even if that value
/// is 0) that quickly can become tedious.
///
/// ```
/// use astroxide::units::*;
/// let long_init = Dimension {length: 1, time: -1, ..Dimension::ZERO};
/// let easy = dim!(length: 1, time: -1);
/// assert_eq!(long_init, easy);
///
/// ```
#[macro_export]
macro_rules! dim {
    ($($field:ident: $value:expr),* $(,)?) => {
        Dimension {
            $($field: $value,)*
            ..Dimension::ZERO
        }
    };
}
pub use dim;

macro_rules! create_base_unit {
    ($name: ident, $symbol: expr, $dimension: expr, $conversion_factor: expr) => {
        #[doc = "Base Unit representation for "]
        #[doc=$symbol]
        pub static $name: BaseUnit = BaseUnit {
            base_dimension: $dimension,
            symbol: $symbol,
            conversion_factor: $conversion_factor,
        };
    };
}

macro_rules! unit_generator {
    ($macro_name: ident, $dimensions: expr) => {
        macro_rules! $macro_name {
            ($name: ident, $symbol: expr, $conversion_factor: expr) => {
                create_base_unit!($name, $symbol, $dimensions, $conversion_factor);
            };
        }
    };
}
unit_generator!(length, dim!(length: 1));
unit_generator!(mass, dim!(mass: 1));
unit_generator!(time, dim!(time: 1));
unit_generator!(temperature, dim!(temperature: 1));
unit_generator!(current, dim!(current: 1));
unit_generator!(angular_distance, dim!(angular_distance: 1));
unit_generator!(solid_angle, dim!(solid_angle: 1));
unit_generator!(luminous_intensity, dim!(luminous_intensity: 1));
unit_generator!(amount_of_substance, dim!(amount_of_substance: 1));
unit_generator!(frequency, dim!(time: -1));
unit_generator!(force, dim!(mass: 1, length: 1, time: -2));
unit_generator!(pressure, dim!(mass: 1, length: -1, time: -2));
unit_generator!(energy, dim!(mass: 1, length: 2, time: -2));
unit_generator!(power, dim!(mass: 1, length: 2, time: -3));
unit_generator!(charge, dim!(time: 1, current: 1));
unit_generator!(voltage, dim!(mass: 1, length: 2, time: -3, current: -1));
unit_generator!(resistance, dim!(mass: 1, length: 2, time: -3, current: -2));
unit_generator!(conductance, dim!(mass: -1, length: -2, time: 3, current: 2));
unit_generator!(capacitance, dim!(mass: -1, length: -2, time: 4, current: 2));
unit_generator!(inductance, dim!(mass: 1, length: 2, time: -2, current: -2));
unit_generator!(magnetic_flux_density, dim!(mass: 1, time: -2, current: -1));
unit_generator!(
    magnetic_flux,
    dim!(mass: 1, length: 2, time: -2, current: -1)
);
unit_generator!(luminous_flux, dim!(luminous_intensity: 1, solid_angle: 1));
unit_generator!(illuminance, dim!(luminous_intensity: 1, length: -2));
unit_generator!(radioactivity, dim!(time: -1));
unit_generator!(absorbed_dose, dim!(length: 2, time: -2));
unit_generator!(equivalent_dose, dim!(length: 2, time: -2));
unit_generator!(catalytic_activity, dim!(time: -1, amount_of_substance: 1));

macro_rules! si {
    ($base_unit: ident, $base_symbol: expr, $base_conversion: expr, $create_macro: ident) => {
        paste! {
            $create_macro!([<QUETTA $base_unit>], concat!("Q", $base_symbol), 1e30 * $base_conversion);
            $create_macro!([<RONNA $base_unit>], concat!("R", $base_symbol), 1e27 * $base_conversion);
            $create_macro!([<YOTTA $base_unit>], concat!("Y", $base_symbol), 1e24 * $base_conversion);
            $create_macro!([<ZETTA $base_unit>], concat!("Z", $base_symbol), 1e21 * $base_conversion);
            $create_macro!([<EXA $base_unit>], concat!("E", $base_symbol), 1e18 * $base_conversion);
            $create_macro!([<PETA $base_unit>], concat!("P", $base_symbol), 1e15 * $base_conversion);
            $create_macro!([<TERA $base_unit>], concat!("T", $base_symbol), 1e12 * $base_conversion);
            $create_macro!([<GIGA $base_unit>], concat!("G", $base_symbol), 1e9 * $base_conversion);
            $create_macro!([<MEGA $base_unit>], concat!("M", $base_symbol), 1e6 * $base_conversion);
            $create_macro!([<KILO $base_unit>], concat!("k", $base_symbol), 1e3 * $base_conversion);
            $create_macro!([<HECTO $base_unit>], concat!("h", $base_symbol), 1e2 * $base_conversion);
            $create_macro!([<DECA $base_unit>], concat!("da", $base_symbol), 1e1 * $base_conversion);
            $create_macro!($base_unit, $base_symbol, $base_conversion);
            $create_macro!([<DECI $base_unit>], concat!("d", $base_symbol), 1e-1 * $base_conversion);
            $create_macro!([<CENTI $base_unit>], concat!("c", $base_symbol), 1e-2 * $base_conversion);
            $create_macro!([<MILLI $base_unit>], concat!("m", $base_symbol), 1e-3 * $base_conversion);
            $create_macro!([<MICRO $base_unit>], concat!("μ", $base_symbol), 1e-6 * $base_conversion);
            $create_macro!([<NANO $base_unit>], concat!("n", $base_symbol), 1e-9 * $base_conversion);
            $create_macro!([<PICO $base_unit>], concat!("p", $base_symbol), 1e-12 * $base_conversion);
            $create_macro!([<FEMTO $base_unit>], concat!("f", $base_symbol), 1e-15 * $base_conversion);
            $create_macro!([<ATTO $base_unit>], concat!("a", $base_symbol), 1e-18 * $base_conversion);
            $create_macro!([<ZEPTO $base_unit>], concat!("z", $base_symbol), 1e-21 * $base_conversion);
            $create_macro!([<YOCTO $base_unit>], concat!("y", $base_symbol), 1e-24 * $base_conversion);
            $create_macro!([<RONTO $base_unit>], concat!("r", $base_symbol), 1e-27 * $base_conversion);
            $create_macro!([<QUECTO $base_unit>], concat!("q", $base_symbol), 1e-30 * $base_conversion);
        }
    };
}

si!(METER, "m", 1., length);
length!(ANGSTROM, "Å", 1e-10); // Adding missing Angstrom
// Astronomical Length Units
si!(ASTRONOMICAL_UNIT, "AU", 1.496e11, length);
si!(LIGHTYEAR, "lyr", 9.5e15, length);
si!(PARSEC, "pc", 3.09e16, length);

// Imperial Length Units
length!(TWIP, "twip", 0.000_017_638_888_888);
length!(THOU, "th", 0.000_025_4);
length!(BARLEYCORN, "barleycorn", 0.008_466_666_666);
length!(INCH, "in", 0.0254);
length!(HAND, "hh", 0.1016);
length!(FOOT, "ft", 0.3048);
length!(YARD, "yd", 0.9144);
length!(CHAIN, "ch", 20.1168);
length!(FURLONG, "fur", 201.168);
length!(MILE, "mi", 1609.344);
length!(LEAGUE, "lea", 4828.032);
length!(FATHOM, "ftm", 1.8288);
length!(CABLE, "cable", 185.2);
length!(NAUTICAL_MILE, "nmi", 1852.);
length!(LINK, "link", 0.201_168);
length!(ROD, "rod", 5.0292);

// Mass
si!(GRAM, "g", 1., mass);

// Astronomical Mass Units
mass!(SOLAR_MASS, "M☉", 1.988_475e33);

// Time
si!(SECOND, "s", 1., time);
si!(MINUTE, "min", 60., time);
si!(HOUR, "hr", 3600., time);

// Metric Temperature units
si!(KELVIN, "K", 1., temperature);

// Angular Distance Units
si!(RADIAN, "rad", 1., angular_distance);
si!(DEGREE, "deg", std::f64::consts::PI / 180., angular_distance);
si!(
    ARCMINUTE,
    "arcmin",
    std::f64::consts::PI / (60. * 180.),
    angular_distance
);
si!(
    ARCSECOND,
    "arcsec",
    std::f64::consts::PI / (3600. * 180.),
    angular_distance
);

// Current Units
si!(AMPERE, "A", 1., current);

// Solid Angle
si!(STERADIAN, "sr", 1., solid_angle);

// Luminous Intensity
si!(CANDELA, "cd", 1., luminous_intensity);

// Amount of substance.
si!(MOL, "mol", 1., amount_of_substance);

// Derived Units
si!(HERTZ, "Hz", 1., frequency);
si!(NEWTON, "N", 1., force);
si!(PASCAL, "Pa", 1., pressure);
si!(JOULE, "J", 1., energy);
si!(ERG, "erg", 1e-7, energy);
si!(WATT, "W", 1., power);
power!(SOLAR_LUM, "L☉", 3.828e26);
si!(COULOMB, "C", 1., charge);
si!(VOLT, "V", 1., voltage);
si!(OHM, "Ω", 1., resistance);
si!(SIEMAN, "S", 1., conductance);
si!(FARAD, "F", 1., capacitance);
si!(HENRY, "H", 1., inductance);
si!(TESLA, "T", 1., magnetic_flux_density);
si!(WEBER, "Wb", 1., magnetic_flux);
si!(LUMEN, "lm", 1., luminous_flux);
si!(LUX, "lx", 1., illuminance);
si!(BECQUEREL, "Bq", 1., radioactivity);
si!(GRAY, "Gy", 1., absorbed_dose);
si!(SIEVERT, "Sv", 1., equivalent_dose);
si!(KATAL, "kat", 1., catalytic_activity);

#[cfg(test)]
mod tests {

    use super::*;
    #[test]
    fn test_unit_factor() {
        let a = 5. * KILOMETER;
        let result = a.unit.calculate_conversion_factor();
        assert!((result - 1000.).abs() < f64::EPSILON);
    }
    #[test]
    fn test_unit_equality() {
        //test that order doens't matter and that unts are actually equal
        let meter_squared = METER * METER;
        let meter_kilometer = METER * KILOMETER;
        let kilometer_squared = KILOMETER * KILOMETER;
        let another_km_squared = KILOMETER * KILOMETER;
        let another_meter_squared = METER * METER;
        assert_eq!(meter_squared, another_meter_squared);
        assert_eq!(kilometer_squared, another_km_squared);
        assert_ne!(meter_kilometer, meter_squared);
        assert_ne!(meter_squared, kilometer_squared);

        let a = METER * KILOMETER;
        let b = KILOMETER * METER;
        assert_eq!(a, b);
    }
    #[test]
    fn test_multiplying_unit() {
        let a = METER * METER; // base unit
        let b = METER * METER * METER; // unit * base_unit
        let c = METER * KILOMETER * SECOND * KILOMETER;
        let answer_a = Unit {
            base_units: { BTreeMap::from([(METER, 2)]) },
        };
        let answer_b = Unit {
            base_units: { BTreeMap::from([(METER, 3)]) },
        };
        let answer_c = Unit {
            base_units: { BTreeMap::from([(METER, 1), (SECOND, 1), (KILOMETER, 2)]) },
        };
        assert_eq!(a, answer_a);
        assert_eq!(b, answer_b);
        assert_eq!(c, answer_c);
    }

    #[test]
    fn test_dividing_units() {
        let a = METER;
        let b = METER;
        let c = a / b;
        assert_eq!(c.base_units, BTreeMap::new());

        let a = METER;
        let b = SECOND;
        let c = a / b;
        assert_eq!(c.base_units, BTreeMap::from([(METER, 1), (SECOND, -1)]));
    }

    #[test]
    fn test_making_quantities() {
        let distance_a = 5. * METER;
        let distance_b = 2. * METER;
        let distance_c = 3. * KILOMETER;
        let meter_distance = distance_a.clone() + distance_b;
        let another_distance = distance_a + distance_c;

        assert!((meter_distance.value - 7.).abs() < f64::EPSILON);
        assert_eq!(meter_distance.unit, METER.as_unit());
        assert!((another_distance.value - 3005.).abs() < f64::EPSILON);
        assert_eq!(another_distance.unit, METER.as_unit());
    }
    #[test]
    fn zero_quantity() {
        let a = 5. * METER;
        let b = 0. * METER;
        let c = a + b;
        assert!((c.value - 5.).abs() < f64::EPSILON);
        assert_eq!(c.unit, METER.as_unit());
    }
    #[test]
    fn adding_and_subtracting_equivalent_units() {
        let a = 5. * METER * 3. * SECOND;
        let b = (3. * HOUR) * (5. * METER);
        let c = a.clone() + b.clone();
        let d = b.clone() - a.clone();
        assert!((c.value - 54015.).abs() < 1e-12);
        assert_eq!(c.unit, a.unit);
        assert!((d.value - 14.995_833_33).abs() < 1e-7);
        assert_eq!(d.unit, b.unit);
    }

    #[test]
    fn test_adding_and_converting_units() {
        let distance_a = 5. * METER;
        let distance_b = 2. * METER;
        let meter_distance = distance_a + distance_b;

        let km_distance = meter_distance.to(KILOMETER);
        let distance = 5. * METER + 1. * KILOMETER;
        assert!((distance.value - 1005.).abs() < f64::EPSILON);
        assert!((meter_distance.value - 7.).abs() < f64::EPSILON);
        assert!((km_distance.value - 0.007).abs() < f64::EPSILON);
    }

    #[test]
    #[should_panic(expected = "Cannot add or subtract")]
    fn test_adding_non_equivalent_units_fails() {
        let distance = 4. * METER;
        let time = 3. * HOUR;
        let _ = distance + time;
    }

    #[test]
    #[should_panic(expected = "Cannot add or subtract")]
    fn test_subtracting_non_equivalent_units_fails() {
        let distance = 4. * METER;
        let time = 3. * HOUR;
        let _ = distance - time;
    }

    #[test]
    fn test_subtracting_units() {
        let distance_a = 5. * KILOMETER;
        let distance_b = 500. * METER;
        let distance_c = CENTIMETER * 500.;
        let test_distance = distance_a - distance_b - distance_c;
        let test_distance_meters = test_distance.to(METER);
        assert!((test_distance.value - 4.495).abs() < f64::EPSILON);
        assert!((test_distance_meters.value - 4495.).abs() < f64::EPSILON);
    }

    #[test]
    fn test_printing_derived_units() {
        let unit = METER * SECOND * KILOMETER * METER;
        let answer = String::from("m² s km");
        let result = unit.to_string();
        let mut answer_strings: Vec<&str> = answer.split(' ').collect();
        let mut result_strings: Vec<&str> = result.split(' ').collect();
        answer_strings.sort_unstable();
        result_strings.sort_unstable();
        assert_eq!(answer_strings, result_strings);
    }

    #[test]
    fn test_printing_quantities() {
        let quant = 5. * MEGAPARSEC / SECOND;
        assert_eq!(quant.to_string(), "5 s⁻¹ Mpc");
    }

    #[test]
    fn test_derived_units() {
        let x = 5. * KILOMETER;
        let y = 10. * SECOND;
        let velocity = x / y;
        let velocity_mh = velocity.to(METER / HOUR);
        assert!((velocity.value - 0.5).abs() < f64::EPSILON);
        assert!((velocity_mh.value - 1_800_000.).abs() < f64::EPSILON);
    }

    #[test]
    fn test_astronomy() {
        let volume = 3. * (MEGAPARSEC * MEGAPARSEC * MEGAPARSEC);
        let volume_gpc3 = volume.to(GIGAPARSEC * GIGAPARSEC * GIGAPARSEC);
        assert!((volume_gpc3.value - 3e-9).abs() < f64::EPSILON);

        let speed = 50. * (METER / SECOND);
        let speed_kmh = speed.to(KILOMETER / HOUR);
        dbg!(&speed_kmh);
        assert!((speed_kmh.value - 180.).abs() < 1e-12);

        let hubble_constant = 70. * (KILOMETER / SECOND) / MEGAPARSEC;
        let weird_hubble = hubble_constant.to(METER / (KILOMETER * HOUR));
        assert!((weird_hubble.value - 8.166_763_81e-12).abs() < 1e-12); // comparing to astropyj
    }
    #[test]
    #[should_panic(expected = "Cannot convert")]
    fn converting_non_equivalent_units_fail() {
        let volume = 5. * METER * METER * METER;
        let _ = volume.to(METER);
    }

    #[test]
    fn testing_angular_conversions() {
        let a = 5. * RADIAN;
        let b = 2. * DEGREE;
        assert!((a.to(DEGREE).value - 286.478_897_565_411_6).abs() < f64::EPSILON);
        assert!((a.to(ARCSECOND).value - 1_031_324.031_235_481_8).abs() < f64::EPSILON);
        assert!((a.to(ARCMINUTE).value - 17_188.733_853_924_7).abs() < f64::EPSILON);

        assert!((b.to(RADIAN).value - 0.034_906_585_039_886_59).abs() < f64::EPSILON);
        assert!((b.to(ARCSECOND).value - 3600. * 2.).abs() < f64::EPSILON);
        assert!((b.to(ARCMINUTE).value - 60. * 2.).abs() < f64::EPSILON);
    }

    #[test]
    fn test_derived_equivalence() {
        let a = 5. * HERTZ;
        let b = a.to((1. / HOUR).unit);
        assert!((b.value - 18000.).abs() < f64::EPSILON);

        let a = 2. * SOLAR_LUM;
        let b = a.to(WATT);
        assert!((b.value - 7.656e26).abs() < f64::EPSILON);
    }

    #[test]
    fn test_factoring_out_h() {
        let a = 1. * MEGAPARSEC;
        let b = a.factor_out_h(0.7, -1);
        let c = a.factor_out_h(0.7, -2);
        assert!((b.cosmo_value.value - 0.7).abs() < f64::EPSILON);
        assert_eq!(b.unit, a.unit);
        assert!((c.cosmo_value.value - 0.49).abs() < f64::EPSILON);
    }

    #[test]
    fn test_assuming_h() {
        let a = 1. * MEGAPARSEC;
        let b = a.factor_out_h(0.7, -1);
        let c = b.factor_in_h(0.7);
        assert!((c.value - a.value).abs() < f64::EPSILON);
        assert_eq!(c.unit, a.unit);
    }

    #[test]
    fn test_adding_subtracting_cosmo_units() {
        let a = CosmoQuantity::new(1., -1, MEGAPARSEC);
        let b = CosmoQuantity::new(0.7, -1, MEGAPARSEC);
        let c = CosmoQuantity::new(2.7, -1, MEGAPARSEC);
        let mul_add = a.clone() + b.clone() + c.clone();
        let answer = CosmoQuantity::new(4.4, -1, MEGAPARSEC);
        assert_eq!(mul_add, answer);

        let mul_sub = a.clone() - b - c;
        let answer = CosmoQuantity::new(-2.4, -1, MEGAPARSEC);
        assert!(mul_sub.approx_eq(&answer, 9));

        let d = CosmoQuantity::new(2., -1, GIGAPARSEC);
        let add_different_units = d + a;
        let answer = CosmoQuantity::new(2.001, -1, GIGAPARSEC);
        assert_eq!(add_different_units, answer);
    }

    #[test]
    #[should_panic(expected = "Cannot add or subtract")]
    fn test_adding_wrong_h_dependencies() {
        let _ = CosmoQuantity::new(1., -1, METER) + CosmoQuantity::new(1., -2, METER);
    }
    #[test]
    #[should_panic(expected = "Cannot add or subtract")]
    fn test_subtracting_wrong_h_dependencies() {
        let _ = CosmoQuantity::new(1., -1, METER) - CosmoQuantity::new(1., -2, METER);
    }

    #[test]
    fn test_multipying_dividing_cosmo_units() {
        let a = CosmoValue::new(1., -1) * MEGAPARSEC;
        let b = CosmoQuantity::new(0.7, -1, MEGAPARSEC);
        let c = CosmoQuantity::new(2.7, -1, MEGAPARSEC);
        let mul_mul = a.clone() * b.clone() * c.clone();
        let answer = CosmoQuantity::new(1.89, -3, MEGAPARSEC * MEGAPARSEC * MEGAPARSEC);
        assert_eq!(mul_mul, answer);

        let mul_div = (a.clone() / b) / c;
        let answer = CosmoQuantity::new(0.529_100_529_1, 1, (1. / MEGAPARSEC).unit);
        assert!(mul_div.approx_eq(&answer, 9));

        let d = CosmoQuantity::new(2., -1, SECOND);
        let different_units = a / d;
        let answer = CosmoQuantity::new(0.5, 0, MEGAPARSEC / SECOND);
        assert_eq!(different_units, answer);
    }

    #[test]
    fn test_converting_units() {
        let a = 1. * MEGAPARSEC;
        let b = a.switch_cosmologies(0.7, 0.6, -1);
        let answer = 1.166_666_666_67 * MEGAPARSEC;
        assert!(b.approx_eq(&answer, 9));
    }

    #[test]
    fn test_print() {
        let little_h = 0.7;
        let plain = 1.78 * MEGAPARSEC;
        let cosmo = plain.factor_out_h(little_h, -1);

        println!("Value assuming h={little_h}: {plain}");
        println!("Value with {}: {}\n", "h".italic(), cosmo);

        let len_1 = 5. * h(-1) * MEGAPARSEC;
        let len_2 = 5. * h(-1) * MEGAPARSEC;
        let len_3 = 5. * h(-1) * MEGAPARSEC;
        let volume = len_1 * len_2 * len_3;
        let volume_gpc3 = volume.to(GIGAPARSEC * GIGAPARSEC * GIGAPARSEC);

        println!("Volume with h: {volume}");
        println!("Volume with h but in GPC: {volume_gpc3}");
        println!(
            "Volume assuming h={}: {}\n",
            little_h,
            volume.factor_in_h(little_h)
        );

        let mass = 5. * SOLAR_MASS;
        let h_mass = mass.factor_out_h(little_h, -1);
        let mass_067 = h_mass.factor_in_h(0.67);

        println!("Mass assuming h=0.7: {mass}");
        println!("Mass with h factored out: {h_mass}");
        println!("Mass assuming h=0.67: {mass_067}");
    }
}
