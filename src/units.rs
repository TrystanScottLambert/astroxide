use fmtastic::Superscript;
use paste::paste;
use std::{
    fmt::{self, Display},
    ops::{Add, Div, Mul, Sub},
};

pub type Result<T> = core::result::Result<T, UnitError>;

#[derive(Debug, Clone)]
pub enum UnitError {
    DifferentDimensions(Vec<Dimension>, Vec<Dimension>),
}

fn are_dimensions_equal(dimensions_a: &[Dimension], dimensions_b: &[Dimension]) -> bool {
    for dim in dimensions_a {
        if !dimensions_b.contains(dim) {
            return false;
        }
    }
    for dim in dimensions_b {
        if !dimensions_a.contains(dim) {
            return false;
        }
    }
    true
}

// Helper function to iterate over the implemented units and remove any that
// have zero in the exponents. If there is nothing left then it returns Unitless
fn clean_up_zero_exponents(units: Vec<ImplBaseUnit>) -> Vec<ImplBaseUnit> {
    let mut new_vec = Vec::new();
    for iu in units {
        if iu.exponent != 0 && iu.base_unit != UNITLESS {
            new_vec.push(iu);
        }
    }
    if new_vec.is_empty() {
        vec![ImplBaseUnit {
            base_unit: UNITLESS,
            exponent: 0,
        }]
    } else {
        new_vec
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BaseDimension {
    Length,
    Mass,
    Time,
    Temperature,
    Current,
    AgularDistance,
    SolidAngle,
    LuminousIntensity,
    AmountOfSubstance,
    Unitless,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dimension {
    pub base: BaseDimension,
    pub exponent: i32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaseUnit {
    pub base_dimension: Dimension,
    pub symbol: &'static str,
    pub conversion_factor: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImplBaseUnit {
    pub base_unit: BaseUnit,
    pub exponent: i32,
}

#[derive(Debug, Clone)]
pub struct Unit {
    pub base_units: Vec<ImplBaseUnit>,
}

#[derive(Debug, Clone)]
pub struct Quantity {
    pub unit: Unit,
    pub value: f64,
}

pub trait UnitLike {
    fn as_unit(&self) -> Unit;
    fn calculate_conversion_factor(&self) -> f64 {
        self.as_unit()
            .base_units
            .iter()
            .map(|iu| iu.base_unit.conversion_factor.powi(iu.exponent))
            .fold(1., |acc, x| acc * x)
    }
}

impl Display for ImplBaseUnit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let repr = match self.exponent {
            0 => String::new(),
            1 => self.base_unit.symbol.to_string(),
            _ => format!("{}{}", self.base_unit.symbol, Superscript(self.exponent)),
        };
        write!(f, "{}", repr)
    }
}

impl Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let string_reprs: Vec<String> =
            self.base_units.iter().map(|iu| format!("{}", iu)).collect();
        write!(f, "{}", string_reprs.join(" "))
    }
}
impl Unit {
    pub fn get_units_list(&self) -> Vec<BaseUnit> {
        self.base_units.iter().map(|iu| iu.base_unit).collect()
    }

    pub fn dimensions(self) -> Vec<Dimension> {
        self.base_units
            .iter()
            .map(|iu| Dimension {
                base: iu.base_unit.base_dimension.base,
                exponent: iu.base_unit.base_dimension.exponent * iu.exponent,
            })
            .collect()
    }
    pub fn invert(&self) -> Self {
        let base_units = self
            .base_units
            .iter()
            .map(|iu| ImplBaseUnit {
                base_unit: iu.base_unit,
                exponent: -iu.exponent,
            })
            .collect();
        Unit { base_units }
    }
}

impl Quantity {
    pub fn to(self, target_unit: impl UnitLike) -> Result<Quantity> {
        if !are_dimensions_equal(
            &self.unit.clone().dimensions(),
            &target_unit.as_unit().dimensions(),
        ) {
            return Err(UnitError::DifferentDimensions(
                self.unit.dimensions(),
                target_unit.as_unit().dimensions(),
            ));
        }
        Ok(Quantity {
            unit: target_unit.as_unit(),
            value: self.value
                * (self.unit.calculate_conversion_factor()
                    / target_unit.calculate_conversion_factor()),
        })
    }
}

impl Display for Quantity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.unit)
    }
}

impl UnitLike for Unit {
    fn as_unit(&self) -> Unit {
        self.clone()
    }
}

impl UnitLike for BaseUnit {
    fn as_unit(&self) -> Unit {
        Unit {
            base_units: vec![ImplBaseUnit {
                base_unit: *self,
                exponent: 1,
            }],
        }
    }
}

impl PartialEq for Unit {
    fn eq(&self, other: &Self) -> bool {
        for iu in self.base_units.clone() {
            if !other.base_units.contains(&iu) {
                return false;
            }
        }
        for iu in other.base_units.clone() {
            if !self.base_units.contains(&iu) {
                return false;
            }
        }
        true
    }
}

impl Div for Dimension {
    type Output = Vec<Dimension>;
    fn div(self, rhs: Self) -> Self::Output {
        if self.base == rhs.base {
            if self.exponent == rhs.exponent {
                vec![Dimension {
                    base: BaseDimension::Unitless,
                    exponent: 0,
                }]
            } else {
                vec![Dimension {
                    base: self.base,
                    exponent: self.exponent - rhs.exponent,
                }]
            }
        } else {
            let negative_rhs = Dimension {
                base: rhs.base,
                exponent: -rhs.exponent,
            };
            vec![self, negative_rhs]
        }
    }
}

impl Div<BaseUnit> for BaseUnit {
    type Output = Unit;
    fn div(self, rhs: BaseUnit) -> Self::Output {
        self.as_unit() / rhs.as_unit()
    }
}

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

        Quantity { value, unit }
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
            value: rhs / self.value,
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
        rhs / self
    }
}

impl Mul for Dimension {
    type Output = Vec<Dimension>;
    fn mul(self, rhs: Self) -> Self::Output {
        if self.base == rhs.base {
            if self.exponent == -rhs.exponent {
                vec![Dimension {
                    base: BaseDimension::Unitless,
                    exponent: 0,
                }]
            } else {
                vec![Dimension {
                    base: self.base,
                    exponent: self.exponent + rhs.exponent,
                }]
            }
        } else {
            vec![self, rhs]
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
        let mut new_impl_units = Vec::new();
        let self_base_units = self.get_units_list();
        let rhs_base_unit_list = rhs.get_units_list();

        for (i, base_unit) in self_base_units.iter().enumerate() {
            if !rhs_base_unit_list.contains(base_unit) {
                new_impl_units.push(*self.base_units.get(i).unwrap())
            }
        }

        for impl_unit in rhs.base_units {
            if self_base_units.contains(&impl_unit.base_unit) {
                for iu in self.base_units.clone() {
                    if iu.base_unit == impl_unit.base_unit {
                        new_impl_units.push(ImplBaseUnit {
                            base_unit: iu.base_unit,
                            exponent: iu.exponent + impl_unit.exponent,
                        })
                    }
                }
            } else {
                new_impl_units.push(impl_unit);
            }
        }
        Unit {
            base_units: clean_up_zero_exponents(new_impl_units),
        }
    }
}

impl Mul for Quantity {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let unit = self.unit * rhs.unit;
        let value = self.value * rhs.value;
        Quantity { value, unit }
    }
}

impl Mul<f64> for BaseUnit {
    type Output = Quantity;
    fn mul(self, rhs: f64) -> Self::Output {
        let implemented_unit = ImplBaseUnit {
            base_unit: self,
            exponent: 1,
        };
        let unit: Unit = Unit {
            base_units: vec![implemented_unit],
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
    type Output = Result<Self>;
    fn add(self, rhs: Self) -> Self::Output {
        let self_dimensions = self.unit.clone().dimensions();
        let rhs_dimensions = rhs.unit.clone().dimensions();
        if !are_dimensions_equal(&self_dimensions, &rhs_dimensions) {
            return Err(UnitError::DifferentDimensions(
                self_dimensions,
                rhs_dimensions,
            ));
        }
        let new_value = (self.value * self.unit.calculate_conversion_factor()
            + rhs.value * rhs.unit.calculate_conversion_factor())
            / self.unit.calculate_conversion_factor();
        Ok(Quantity {
            unit: self.unit.clone(),
            value: new_value,
        })
    }
}
impl Add<Quantity> for Result<Quantity> {
    type Output = Self;
    fn add(self, rhs: Quantity) -> Self::Output {
        let x = self?;
        x + rhs
    }
}

impl Add<Result<Quantity>> for Quantity {
    type Output = Result<Quantity>;
    fn add(self, rhs: Result<Quantity>) -> Self::Output {
        rhs + self
    }
}

impl Sub for Quantity {
    type Output = Result<Self>;
    fn sub(self, rhs: Self) -> Self::Output {
        let self_dimensions = self.unit.clone().dimensions();
        let rhs_dimensions = rhs.unit.clone().dimensions();
        if !are_dimensions_equal(&self_dimensions, &rhs_dimensions) {
            return Err(UnitError::DifferentDimensions(
                self_dimensions,
                rhs_dimensions,
            ));
        }
        let new_value = (self.value * self.unit.calculate_conversion_factor()
            - rhs.value * rhs.unit.calculate_conversion_factor())
            / self.unit.calculate_conversion_factor();
        Ok(Quantity {
            unit: self.unit.clone(),
            value: new_value,
        })
    }
}

impl Sub<Quantity> for Result<Quantity> {
    type Output = Self;
    fn sub(self, rhs: Quantity) -> Self::Output {
        let x = self?;
        x - rhs
    }
}

impl Sub<Result<Quantity>> for Quantity {
    type Output = Result<Quantity>;
    fn sub(self, rhs: Result<Quantity>) -> Self::Output {
        let x = rhs?;
        self - x
    }
}

macro_rules! create_base_unit {
    ($name: ident, $symbol: expr, $dimension: expr, $conversion_factor: expr, $exponent: expr) => {
        pub static $name: BaseUnit = BaseUnit {
            base_dimension: Dimension {
                base: $dimension,
                exponent: $exponent,
            },
            symbol: $symbol,
            conversion_factor: $conversion_factor,
        };
    };
}

create_base_unit!(UNITLESS, "", BaseDimension::Unitless, 0., 0);

macro_rules! create_length_unit {
    ($name: ident, $symbol: expr, $conversion_factor: expr) => {
        create_base_unit!($name, $symbol, BaseDimension::Length, $conversion_factor, 1);
    };
}
macro_rules! create_mass_unit {
    ($name: ident, $symbol: expr, $conversion_factor: expr) => {
        create_base_unit!($name, $symbol, BaseDimension::Mass, $conversion_factor, 1);
    };
}
macro_rules! create_time_unit {
    ($name: ident, $symbol: expr, $conversion_factor: expr) => {
        create_base_unit!($name, $symbol, BaseDimension::Time, $conversion_factor, 1);
    };
}

macro_rules! create_temperature_unit {
    ($name: ident, $symbol: expr, $conversion_factor: expr) => {
        create_base_unit!(
            $name,
            $symbol,
            BaseDimension::Temperature,
            $conversion_factor,
            1
        );
    };
}

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
            $create_macro!([<MICRO $base_unit>], concat!("u", $base_symbol), 1e-6 * $base_conversion);
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

si!(METER, "m", 1., create_length_unit);
// Astronomical Length Units
si!(ASTRONOMICAL_UNIT, "AU", 1.496e11, create_length_unit);
si!(LIGHTYEAR, "lyr", 9.5e15, create_length_unit);
si!(PARSEC, "pc", 3.09e16, create_length_unit);

// Imperial Length Units
create_length_unit!(TWIP, "twip", 0.000017638888888);
create_length_unit!(THOU, "th", 0.0000254);
create_length_unit!(BARLEYCORN, "barelycorn", 0.008466666666);
create_length_unit!(INCH, "in", 0.0254);
create_length_unit!(HAND, "hh", 0.1016);
create_length_unit!(FOOT, "ft", 0.3048);
create_length_unit!(YARD, "yd", 0.9144);
create_length_unit!(CHAIN, "ch", 20.1168);
create_length_unit!(FURLONG, "fur", 201.168);
create_length_unit!(MILE, "mi", 1609.344);
create_length_unit!(LEAGUE, "lea", 4828.032);
create_length_unit!(FATHOM, "ftm", 1.852);
create_length_unit!(CABLE, "cable", 185.2);
create_length_unit!(NAUTICAL_MILE, "nmi", 185.2);
create_length_unit!(LINK, "link", 0.201168);
create_length_unit!(ROD, "rod", 5.0292);

// Mass
si!(GRAM, "g", 1., create_mass_unit);

// Astronomical Mass Units
create_mass_unit!(SOLAR_MASS, "msun", 1.988475e33);

// Time
si!(SECOND, "s", 1., create_time_unit);
si!(MINUTE, "min", 60., create_time_unit);
si!(HOUR, "hr", 3600., create_time_unit);

// Metric Temperature units
si!(KELVIN, "K", 1., create_temperature_unit);

#[cfg(test)]
mod tests {

    use super::*;
    #[test]
    fn test_unit_factor() {
        let a = 5. * KILOMETER;
        let result = a.unit.calculate_conversion_factor();
        assert_eq!(result, 1000.);
    }
    #[test]
    fn test_unit_equality() {
        //test that order doens't matter and that unts are actually equal
        let a = METER * METER;
        let b = METER * KILOMETER;
        let c = KILOMETER * KILOMETER;
        let d = KILOMETER * KILOMETER;
        let e = METER * METER;
        assert_eq!(a, e);
        assert_eq!(c, d);
        assert_ne!(b, a);
        assert_ne!(a, c);

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
            base_units: vec![ImplBaseUnit {
                base_unit: METER,
                exponent: 2,
            }],
        };
        let answer_b = Unit {
            base_units: vec![ImplBaseUnit {
                base_unit: METER,
                exponent: 3,
            }],
        };
        let answer_c = Unit {
            base_units: vec![
                ImplBaseUnit {
                    base_unit: METER,
                    exponent: 1,
                },
                ImplBaseUnit {
                    base_unit: SECOND,
                    exponent: 1,
                },
                ImplBaseUnit {
                    base_unit: KILOMETER,
                    exponent: 2,
                },
            ],
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
        assert_eq!(
            c.base_units,
            vec![ImplBaseUnit {
                base_unit: UNITLESS,
                exponent: 0,
            }]
        );

        let a = METER;
        let b = SECOND;
        let c = a / b;
        let answer_simple = vec![
            ImplBaseUnit {
                base_unit: METER,
                exponent: 1,
            },
            ImplBaseUnit {
                base_unit: SECOND,
                exponent: -1,
            },
        ];
        assert_eq!(c.base_units, answer_simple);
    }

    #[test]
    fn test_making_quantities() {
        let distance_a = 5. * METER;
        let distance_b = 2. * METER;
        let distance_c = 3. * KILOMETER;
        let meter_distance = (distance_a.clone() + distance_b).unwrap();
        let another_distance = (distance_a + distance_c).unwrap();

        assert_eq!(meter_distance.value, 7.);
        assert_eq!(meter_distance.unit, METER.as_unit());
        assert_eq!(another_distance.value, 3005.);
        assert_eq!(another_distance.unit, METER.as_unit());
    }
    #[test]
    fn zero_quantity() {
        let a = 5. * METER;
        let b = 0. * METER;
        let c = (a + b).unwrap();
        assert_eq!(c.value, 5.);
        assert_eq!(c.unit, METER.as_unit());
    }
    #[test]
    fn adding_and_subtracting_equivalent_units() {
        let a = 5. * METER * 3. * SECOND;
        let b = (3. * HOUR) * (5. * METER);
        let c = a.clone() + b.clone();
        let d = b.clone() - a.clone();
        assert!((c.clone().unwrap().value - 54015.).abs() < 1e-12);
        assert_eq!(c.unwrap().unit, a.unit);
        assert!((d.clone().unwrap().value - 14.99583333).abs() < 1e-7);
        assert_eq!(d.unwrap().unit, b.unit);
    }

    #[test]
    fn test_adding_and_converting_units() {
        let distance_a = 5. * METER;
        let distance_b = 2. * METER;
        let meter_distance = (distance_a + distance_b).unwrap();

        let km_distance = meter_distance.clone().to(KILOMETER).unwrap();
        let distance = (5. * METER + 1. * KILOMETER).unwrap();
        assert_eq!(distance.value, 1005.);
        assert_eq!(meter_distance.value, 7.);
        assert_eq!(km_distance.value, 0.007);
    }

    #[test]
    fn test_adding_and_subtracting_non_equivalent_units_fails() {
        let distance = 4. * METER;
        let time = 3. * HOUR;
        let thing = distance + time;
        assert!(thing.is_err())
    }

    #[test]
    fn test_subtracting_units() {
        let distance_a = 5. * KILOMETER;
        let distance_b = 500. * METER;
        let distance_c = CENTIMETER * 500.;
        let test_distance = (distance_a - distance_b - distance_c).unwrap();
        let test_distance_meters = test_distance.clone().to(METER).unwrap();
        assert_eq!(test_distance.value, 4.495);
        assert_eq!(test_distance_meters.value, 4495.);
    }

    #[test]
    fn test_printing_derived_units() {
        let unit = METER * SECOND * KILOMETER * METER;
        let answer = String::from("m² s km");
        let result = unit.to_string();
        let mut answer_strings: Vec<&str> = answer.split(" ").collect();
        let mut result_strings: Vec<&str> = result.split(" ").collect();
        answer_strings.sort();
        result_strings.sort();
        assert_eq!(answer_strings, result_strings)
    }

    #[test]
    fn test_printing_quantities() {
        let quant = 5. * MEGAPARSEC / SECOND;
        assert_eq!(quant.to_string(), "5 Mpc s⁻¹")
    }

    #[test]
    fn test_impl_string_repr() {
        let a = ImplBaseUnit {
            base_unit: MEGAPARSEC,
            exponent: 2,
        };
        assert_eq!(format!("{}", a), String::from("Mpc²"))
    }

    #[test]
    fn test_derived_units() {
        let x = 5. * KILOMETER;
        let y = 10. * SECOND;
        let velocity = x / y;
        let velocity_mh = velocity.clone().to(METER / HOUR).unwrap();
        assert_eq!(velocity.value, 0.5);
        assert_eq!(velocity_mh.value, 1800000.)
    }

    #[test]
    fn test_astronomy() {
        let volume = 3. * (MEGAPARSEC * MEGAPARSEC * MEGAPARSEC);
        let volume_gpc3 = volume.to(GIGAPARSEC * GIGAPARSEC * GIGAPARSEC).unwrap();
        assert!((volume_gpc3.value - 3e-9).abs() < f64::EPSILON);

        let speed = 50. * (METER / SECOND);
        let speed_kmh = speed.to(KILOMETER / HOUR).unwrap();
        dbg!(&speed_kmh);
        assert!((speed_kmh.value - 180.).abs() < 1e-12);

        let hubble_constant = 70. * (KILOMETER / SECOND) / MEGAPARSEC;
        let weird_hubble = hubble_constant.to(METER / (KILOMETER * HOUR)).unwrap();
        assert!((weird_hubble.value - 8.16676381e-12).abs() < 1e-12); // comparing to astropyj
    }
    #[test]
    fn converting_non_equivalent_units_fail() {
        let volume = 5. * METER * METER * METER;
        let x = volume.to(METER);
        assert!(x.is_err());
    }
    #[test]
    fn dimensions_are_equal_function() {
        let dim1 = Dimension {
            base: BaseDimension::Length,
            exponent: 1,
        };
        let dim2 = Dimension {
            base: BaseDimension::Length,
            exponent: 2,
        };
        let dim3 = Dimension {
            base: BaseDimension::Mass,
            exponent: 1,
        };
        let dim4 = Dimension {
            base: BaseDimension::Time,
            exponent: 1,
        };
        let a = vec![dim1, dim3];
        let b = vec![dim3, dim1];
        let c = vec![dim2, dim3];
        assert!(are_dimensions_equal(&a, &b));
        assert!(!are_dimensions_equal(&a, &c));
        assert!(are_dimensions_equal(&[dim1], &[dim1]));
        assert!(are_dimensions_equal(&[dim1, dim4], &[dim4, dim1]));
    }
}
