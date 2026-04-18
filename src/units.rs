use colored::Colorize;
use fmtastic::Superscript;
use paste::paste;
use std::cmp::Ordering;
use std::{
    collections::BTreeMap,
    fmt::{self, Display},
    ops::{Add, Div, Mul, Neg, Sub},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Dimension {
    pub length: i32,
    pub mass: i32,
    pub time: i32,
    pub temperature: i32,
    pub current: i32,
    pub angular_distance: i32,
    pub solid_angle: i32,
    pub luminous_intensity: i32,
    pub amount_of_substance: i32,
}

impl Dimension {
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

#[derive(Debug, Clone, Copy)]
pub struct BaseUnit {
    pub base_dimension: Dimension,
    pub symbol: &'static str,
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

#[derive(Debug, Clone)]
pub struct Unit {
    pub base_units: BTreeMap<BaseUnit, i32>,
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
    pub fn dimensions(&self) -> Dimension {
        let mut full_dimension = Dimension::new();
        for (base_unit, &exponent) in &self.base_units {
            full_dimension = full_dimension + (base_unit.base_dimension * exponent)
        }
        full_dimension
    }
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
    pub fn new(value: f64, unit: impl UnitLike) -> Self {
        Self {
            value,
            unit: unit.as_unit(),
        }
    }
    pub fn to(&self, target_unit: impl UnitLike) -> Self {
        if self.unit.dimensions() != target_unit.as_unit().dimensions() {
            panic!(
                "Cannot convert {} to {} since they have different dimensions.",
                self.unit,
                target_unit.as_unit()
            )
        }

        Quantity {
            unit: target_unit.as_unit(),
            value: self.value
                * (self.unit.calculate_conversion_factor()
                    / target_unit.calculate_conversion_factor()),
        }
    }
    pub fn factor_out_h(&self, h_value: f64, h_dependency: i32) -> CosmoQuantity {
        CosmoQuantity::new(
            self.value / (h_value.powi(h_dependency)),
            h_dependency,
            self.unit.clone(),
        )
    }
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
        Quantity { value, unit }
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
        if self.unit.dimensions() != rhs.unit.dimensions() {
            panic!(
                "Cannot add units {} and {} since they have different dimensions.",
                self.unit, rhs.unit,
            )
        }
        let new_value = (self.value * self.unit.calculate_conversion_factor()
            + rhs.value * rhs.unit.calculate_conversion_factor())
            / self.unit.calculate_conversion_factor();
        Quantity {
            unit: self.unit.clone(),
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
        if self.unit.dimensions() != rhs.unit.dimensions() {
            panic!(
                "Cannot subtract units {} and {} since they have different dimensions.",
                self.unit, rhs.unit,
            )
        }
        self + (-rhs)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CosmoValue {
    pub value: f64,
    pub h_dependency: i32,
}

impl CosmoValue {
    pub fn approx_eq(&self, other: &Self) -> bool {
        (self.value - other.value).abs() < 1e-9 && self.h_dependency == other.h_dependency
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

#[allow(non_camel_case_types)]
pub struct h(i32);

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

#[derive(Debug, Clone, PartialEq)]
pub struct CosmoQuantity {
    pub cosmo_value: CosmoValue,
    pub unit: Unit,
}

impl CosmoQuantity {
    pub fn new(value: f64, h_dependency: i32, unit: impl UnitLike) -> Self {
        Self {
            cosmo_value: CosmoValue {
                value,
                h_dependency,
            },
            unit: unit.as_unit(),
        }
    }
    pub fn approx_eq(&self, other: Self) -> bool {
        self.cosmo_value.approx_eq(&other.cosmo_value) && self.unit == other.unit
    }
    pub fn factor_in_h(&self, h_value: f64) -> Quantity {
        Quantity {
            unit: self.unit.clone(),
            value: self.cosmo_value.value * (h_value).powi(self.cosmo_value.h_dependency),
        }
    }

    pub fn ignore_h(&self) -> Quantity {
        Quantity {
            unit: self.unit.clone(),
            value: self.cosmo_value.value,
        }
    }
    pub fn invert(&self) -> CosmoQuantity {
        CosmoQuantity {
            cosmo_value: CosmoValue {
                value: 1. / self.cosmo_value.value,
                h_dependency: -self.cosmo_value.h_dependency,
            },
            unit: self.unit.clone().invert(),
        }
    }
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
            unit: self.unit.clone(),
        }
    }
}

impl Add for CosmoQuantity {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        if self.cosmo_value.h_dependency != rhs.cosmo_value.h_dependency {
            panic!(
                "{} cannot be added to {} because of differing h dependencies.",
                self, rhs
            )
        }
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

macro_rules! dim {
    ($($field:ident: $value:expr),* $(,)?) => {
        Dimension {
            $($field: $value,)*
            ..Dimension::ZERO
        }
    };
}

macro_rules! create_base_unit {
    ($name: ident, $symbol: expr, $dimension: expr, $conversion_factor: expr) => {
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

si!(METER, "m", 1., length);
length!(ANGSTROM, "Å", 1e-10); // Adding missing Angstrom
// Astronomical Length Units
si!(ASTRONOMICAL_UNIT, "AU", 1.496e11, length);
si!(LIGHTYEAR, "lyr", 9.5e15, length);
si!(PARSEC, "pc", 3.09e16, length);

// Imperial Length Units
length!(TWIP, "twip", 0.000017638888888);
length!(THOU, "th", 0.0000254);
length!(BARLEYCORN, "barleycorn", 0.008466666666);
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
length!(LINK, "link", 0.201168);
length!(ROD, "rod", 5.0292);

// Mass
si!(GRAM, "g", 1., mass);

// Astronomical Mass Units
mass!(SOLAR_MASS, "msun", 1.988475e33);

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
power!(SOLAR_LUM, "Lsun", 3.828e26);
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

        assert_eq!(meter_distance.value, 7.);
        assert_eq!(meter_distance.unit, METER.as_unit());
        assert_eq!(another_distance.value, 3005.);
        assert_eq!(another_distance.unit, METER.as_unit());
    }
    #[test]
    fn zero_quantity() {
        let a = 5. * METER;
        let b = 0. * METER;
        let c = a + b;
        assert_eq!(c.value, 5.);
        assert_eq!(c.unit, METER.as_unit());
    }
    #[test]
    fn adding_and_subtracting_equivalent_units() {
        let a = 5. * METER * 3. * SECOND;
        let b = (3. * HOUR) * (5. * METER);
        let c = a.clone() + b.clone();
        let d = b.clone() - a.clone();
        assert!((c.clone().value - 54015.).abs() < 1e-12);
        assert_eq!(c.unit, a.unit);
        assert!((d.clone().value - 14.99583333).abs() < 1e-7);
        assert_eq!(d.unit, b.unit);
    }

    #[test]
    fn test_adding_and_converting_units() {
        let distance_a = 5. * METER;
        let distance_b = 2. * METER;
        let meter_distance = distance_a + distance_b;

        let km_distance = meter_distance.to(KILOMETER);
        let distance = 5. * METER + 1. * KILOMETER;
        assert_eq!(distance.value, 1005.);
        assert_eq!(meter_distance.value, 7.);
        assert_eq!(km_distance.value, 0.007);
    }

    #[test]
    #[should_panic(expected = "Cannot add")]
    fn test_adding_non_equivalent_units_fails() {
        let distance = 4. * METER;
        let time = 3. * HOUR;
        let _ = distance + time;
    }

    #[test]
    #[should_panic(expected = "Cannot subtract")]
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
        assert_eq!(quant.to_string(), "5 s⁻¹ Mpc")
    }

    #[test]
    fn test_derived_units() {
        let x = 5. * KILOMETER;
        let y = 10. * SECOND;
        let velocity = x / y;
        let velocity_mh = velocity.to(METER / HOUR);
        assert_eq!(velocity.value, 0.5);
        assert_eq!(velocity_mh.value, 1800000.)
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
        assert!((weird_hubble.value - 8.16676381e-12).abs() < 1e-12); // comparing to astropyj
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
        assert_eq!(a.to(DEGREE).value, 286.4788975654116);
        assert_eq!(a.to(ARCSECOND).value, 1031324.0312354818);
        assert_eq!(a.to(ARCMINUTE).value, 17188.7338539247);

        assert_eq!(b.to(RADIAN).value, 0.03490658503988659);
        assert_eq!(b.to(ARCSECOND).value, 3600. * 2.);
        assert_eq!(b.to(ARCMINUTE).value, 60. * 2.);
    }

    #[test]
    fn test_derived_equivalence() {
        let a = 5. * HERTZ;
        let b = a.to((1. / HOUR).unit);
        assert_eq!(b.value, 18000.);

        let a = 2. * SOLAR_LUM;
        let b = a.to(WATT);
        assert_eq!(b.value, 7.656e26);
    }

    #[test]
    fn test_factoring_out_h() {
        let a = 1. * MEGAPARSEC;
        let b = a.factor_out_h(0.7, -1);
        let c = a.factor_out_h(0.7, -2);
        assert_eq!(b.cosmo_value.value, 0.7);
        assert_eq!(b.unit, a.unit);
        assert!((c.cosmo_value.value - 0.49).abs() < f64::EPSILON);
    }

    #[test]
    fn test_assuming_h() {
        let a = 1. * MEGAPARSEC;
        let b = a.factor_out_h(0.7, -1);
        let c = b.factor_in_h(0.7);
        assert_eq!(c.value, a.value);
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

        let mul_sub = a.clone() - b.clone() - c.clone();
        let answer = CosmoQuantity::new(-2.4, -1, MEGAPARSEC);
        assert!(mul_sub.approx_eq(answer));

        let d = CosmoQuantity::new(2., -1, GIGAPARSEC);
        let add_different_units = d.clone() + a.clone();
        let answer = CosmoQuantity::new(2.001, -1, GIGAPARSEC);
        assert_eq!(add_different_units, answer);
    }

    #[test]
    #[should_panic]
    fn test_adding_wrong_h_dependencies() {
        let _ = CosmoQuantity::new(1., -1, METER) + CosmoQuantity::new(1., -2, METER);
    }
    #[test]
    #[should_panic]
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

        let mul_div = (a.clone() / b.clone()) / c.clone();
        let answer = CosmoQuantity::new(0.5291005291, 1, (1. / MEGAPARSEC).unit);
        assert!(mul_div.approx_eq(answer));

        let d = CosmoQuantity::new(2., -1, SECOND);
        let different_units = a / d;
        let answer = CosmoQuantity::new(0.5, 0, MEGAPARSEC / SECOND);
        assert_eq!(different_units, answer);
    }

    #[test]
    fn test_converting_units() {
        let a = 0.7 * MEGAPARSEC;
        let b = a.switch_cosmologies(0.7, 0.6, -1);
        let answer = Quantity::new(0.6, MEGAPARSEC);
        assert_eq!(b.value, 0.6)
    }

    #[test]
    // #[should_panic]
    fn test_print() {
        let little_h = 0.7;
        let plain = 1.78 * MEGAPARSEC;
        let cosmo = plain.factor_out_h(little_h, -1);

        println!("Value assuming h={little_h}: {plain}");
        println!("Value with {}: {}", "h".italic(), cosmo);

        let len_1 = 5. * h(-1) * MEGAPARSEC;
        let len_2 = 5. * h(-1) * MEGAPARSEC;
        let len_3 = 5. * h(-1) * MEGAPARSEC;
        let volume = len_1 * len_2 * len_3;
        let volume_gpc3 = volume.to(GIGAPARSEC * GIGAPARSEC * GIGAPARSEC);

        println!("Volume with h: {volume}");
        println!("Volume with h but in GPC: {volume_gpc3}");
        println!(
            "Volume assuming h={}: {}",
            little_h,
            volume.factor_in_h(little_h)
        );

        let mass = 5. * SOLAR_MASS;
        let h_mass = mass.factor_out_h(little_h, -1);
        let mass_067 = h_mass.factor_in_h(0.67);

        println!("Mass assuming h=0.7: {mass}");
        println!("Mass with h factored out: {h_mass}");
        println!("Mass assuming h=0.67: {mass_067}");

        panic!()
    }
}
