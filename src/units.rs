use fmtastic::Superscript;
use std::{
    collections::HashMap,
    ops::{Add, Div, Mul, Sub},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BaseDimension {
    LENGTH,
    MASS,
    TIME,
    TEMPERATURE,
    UNITLESS,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Dimension {
    pub base: BaseDimension,
    pub exponent: i32,
}
impl Mul for Dimension {
    type Output = Vec<Dimension>;
    fn mul(self, rhs: Self) -> Self::Output {
        if self.base == rhs.base {
            if self.exponent == -rhs.exponent {
                vec![Dimension {
                    base: BaseDimension::UNITLESS,
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

impl Div for Dimension {
    type Output = Vec<Dimension>;
    fn div(self, rhs: Self) -> Self::Output {
        if self.base == rhs.base {
            if self.exponent == rhs.exponent {
                vec![Dimension {
                    base: BaseDimension::UNITLESS,
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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaseUnit {
    pub base_dimension: Dimension,
    pub symbol: &'static str,
    pub conversion_factor: f64,
}

#[derive(Debug, Clone)]
pub struct DerivedUnit {
    pub base_units: Vec<BaseUnit>,
    pub symbol: String,
}

impl Mul<BaseUnit> for f64 {
    type Output = BaseQuantity;
    fn mul(self, rhs: BaseUnit) -> Self::Output {
        BaseQuantity {
            unit: rhs,
            value: self,
        }
    }
}

impl Mul<f64> for BaseUnit {
    type Output = BaseQuantity;
    fn mul(self, rhs: f64) -> Self::Output {
        rhs * self
    }
}

pub fn print_unit_from_units(units: Vec<BaseUnit>) -> String {
    let mut unit_count: HashMap<&'static str, i32> = HashMap::new();
    let symbols: Vec<&'static str> = units.iter().map(|u| u.symbol).collect();
    let exponents: Vec<i32> = units.iter().map(|u| u.base_dimension.exponent).collect();

    for (symbol, exp) in symbols.into_iter().zip(exponents) {
        *unit_count.entry(symbol).or_insert(0) += exp;
    }
    let mut output = String::new();
    for (&k, v) in unit_count.iter() {
        let current_unit = format!("{}{}", k, Superscript(*v));
        output.push_str(&current_unit);
    }
    output
}

#[derive(Debug)]
pub struct BaseQuantity {
    pub unit: BaseUnit,
    pub value: f64,
}

impl BaseQuantity {
    pub fn to(&self, unit: BaseUnit) -> BaseQuantity {
        if unit.base_dimension != self.unit.base_dimension {
            panic!("Cannot convert to this unit. Diffent dimensions.")
        } else {
            let converted_value =
                self.value * (self.unit.conversion_factor / unit.conversion_factor);
            BaseQuantity {
                unit,
                value: converted_value,
            }
        }
    }
}

impl Add for BaseQuantity {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        if self.unit.base_dimension == rhs.unit.base_dimension {
            let base_conversion =
                self.value * self.unit.conversion_factor + rhs.value * rhs.unit.conversion_factor;
            let lhs_converted = base_conversion / self.unit.conversion_factor;
            BaseQuantity {
                unit: self.unit,
                value: lhs_converted,
            }
        } else {
            panic!("THESE UNITS CANNOT BE ADDED")
        }
    }
}
impl Sub for BaseQuantity {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        if self.unit.base_dimension == rhs.unit.base_dimension {
            let base_conversion =
                self.value * self.unit.conversion_factor - rhs.value * rhs.unit.conversion_factor;
            let lhs_converted = base_conversion / self.unit.conversion_factor;
            BaseQuantity {
                unit: self.unit,
                value: lhs_converted,
            }
        } else {
            panic!("THESE UNITS CANNOT BE ADDED")
        }
    }
}

macro_rules! create_base_unit {
    ($name: ident, $symbol: expr, $dimension: expr, $conversion_factor: expr, $exponent: expr) => {
        pub static $name: BaseUnit = BaseUnit {
            base_dimension: $dimension,
            symbol: $symbol,
            conversion_factor: $conversion_factor,
            exponent: $exponent,
        };
    };
}

create_base_unit!(UNITLESS, "", BaseDimension::UNITLESS, 0., 0);

macro_rules! create_length_unit {
    ($name: ident, $symbol: expr, $conversion_factor: expr) => {
        create_base_unit!($name, $symbol, BaseDimension::LENGTH, $conversion_factor, 1);
    };
}
macro_rules! create_mass_unit {
    ($name: ident, $symbol: expr, $conversion_factor: expr) => {
        create_base_unit!($name, $symbol, BaseDimension::MASS, $conversion_factor, 1);
    };
}
macro_rules! create_time_unit {
    ($name: ident, $symbol: expr, $conversion_factor: expr) => {
        create_base_unit!($name, $symbol, BaseDimension::TIME, $conversion_factor, 1);
    };
}

macro_rules! create_temperature_unit {
    ($name: ident, $symbol: expr, $conversion_factor: expr) => {
        create_base_unit!(
            $name,
            $symbol,
            BaseDimension::TEMPERATURE,
            $conversion_factor,
            1
        );
    };
}

create_length_unit!(YOTTAMETER, "Ym", 1e24);
create_length_unit!(ZETTAMETER, "Zm", 1e21);
create_length_unit!(EXAMETER, "Em", 1e18);
create_length_unit!(PETAMETER, "Pm", 1e15);
create_length_unit!(TERAMETER, "Tm", 1e12);
create_length_unit!(GIGAMETER, "Gm", 1e9);
create_length_unit!(MEGAMETER, "Mm", 1e6);
create_length_unit!(KILOMETER, "km", 1e3);
create_length_unit!(HECTOMETER, "hm", 1e2);
create_length_unit!(DEKAMETER, "dam", 1e1);
create_length_unit!(METER, "m", 1.);
create_length_unit!(DECIMETER, "dm", 1e-1);
create_length_unit!(CENTIMETER, "cm", 1e-2);
create_length_unit!(MILLIMETER, "mm", 1e-3);
create_length_unit!(MICROMETER, "μm", 1e-6);
create_length_unit!(NANOMETER, "nm", 1e-9);
create_length_unit!(PICOMETER, "pm", 1e-12);
create_length_unit!(FEMTOMETER, "fm", 1e-15);
create_length_unit!(ATTOMETER, "am", 1e-18);
create_length_unit!(ZEPTOMETER, "zm", 1e-21);
create_length_unit!(YOCTOMETER, "ym", 1e-24);

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

// Astronomical Length Units
create_length_unit!(ASTRONOMICAL_UNIT, "AU", 1.496e11);
create_length_unit!(LIGHTYEAR, "lyr", 9.5e15);
create_length_unit!(PARSEC, "pc", 3.09e16);
create_length_unit!(KILO_PARSEC, "kpc", 1e3 * 3.09e16);
create_length_unit!(MEGA_PARSEC, "Mpc", 1e6 * 3.09e16);
create_length_unit!(GIGA_PARSEC, "Gpc", 1e9 * 3.09e16);

// Metric Mass Units
create_mass_unit!(YOTTAGRAM, "Yg", 1e24);
create_mass_unit!(ZETTAGRAM, "Zg", 1e21);
create_mass_unit!(EXAGRAM, "Eg", 1e18);
create_mass_unit!(PETAGRAM, "Pg", 1e15);
create_mass_unit!(TERAGRAM, "Tg", 1e12);
create_mass_unit!(GIGAGRAM, "Gg", 1e9);
create_mass_unit!(MEGAGRAM, "Mg", 1e6);
create_mass_unit!(KILOGRAM, "kg", 1e3);
create_mass_unit!(HECTOGRAM, "hg", 1e2);
create_mass_unit!(DEKAGRAM, "dag", 1e1);
create_mass_unit!(GRAM, "g", 1.);
create_mass_unit!(DECIGRAM, "dg", 1e-1);
create_mass_unit!(CENTIGRAM, "cg", 1e-2);
create_mass_unit!(MILLIGRAM, "mg", 1e-3);
create_mass_unit!(MICROGRAM, "μg", 1e-6);
create_mass_unit!(NANOGRAM, "ng", 1e-9);
create_mass_unit!(PICOGRAM, "pg", 1e-12);
create_mass_unit!(FEMTOGRAM, "fg", 1e-15);
create_mass_unit!(ATTOGRAM, "ag", 1e-18);
create_mass_unit!(ZEPTOGRAM, "zg", 1e-21);
create_mass_unit!(YOCTOGRAM, "yg", 1e-24);

// Astronomical Mass Units
create_mass_unit!(SOLAR_MASS, "msun", 1.988475e33);

// Metric Time Units
create_time_unit!(YOTTASECOND, "Ys", 1e24);
create_time_unit!(ZETTASECOND, "Zs", 1e21);
create_time_unit!(EXASECOND, "Es", 1e18);
create_time_unit!(PETASECOND, "Ps", 1e15);
create_time_unit!(TERASECOND, "Ts", 1e12);
create_time_unit!(GIGASECOND, "Gs", 1e9);
create_time_unit!(MEGASECOND, "Ms", 1e6);
create_time_unit!(KILOSECOND, "ks", 1e3);
create_time_unit!(HECTOSECOND, "hs", 1e2);
create_time_unit!(DEKASECOND, "das", 1e1);
create_time_unit!(SECOND, "s", 1.);
create_time_unit!(DECISECOND, "ds", 1e-1);
create_time_unit!(CENTISECOND, "cs", 1e-2);
create_time_unit!(MILLISECOND, "ms", 1e-3);
create_time_unit!(MICROSECOND, "μs", 1e-6);
create_time_unit!(NANOSECOND, "ns", 1e-9);
create_time_unit!(PICOSECOND, "ps", 1e-12);
create_time_unit!(FEMTOSECOND, "fs", 1e-15);
create_time_unit!(ATTOSECOND, "as", 1e-18);
create_time_unit!(ZEPTOSECOND, "zs", 1e-21);
create_time_unit!(YOCTOSECOND, "ys", 1e-24);

create_time_unit!(MINUTE, "m", 60.);
create_time_unit!(HOUR, "h", 3600.);

// Metric Temperature units
create_temperature_unit!(YOTTAKELVIN, "YK", 1e24);
create_temperature_unit!(ZETTAKELVIN, "ZK", 1e21);
create_temperature_unit!(EXAKELVIN, "EK", 1e18);
create_temperature_unit!(PETAKELVIN, "PK", 1e15);
create_temperature_unit!(TERAKELVIN, "TK", 1e12);
create_temperature_unit!(GIGAKELVIN, "GK", 1e9);
create_temperature_unit!(MEGAKELVIN, "MK", 1e6);
create_temperature_unit!(KILOKELVIN, "kK", 1e3);
create_temperature_unit!(HECTOKELVIN, "hK", 1e2);
create_temperature_unit!(DEKAKELVIN, "daK", 1e1);
create_temperature_unit!(KELVIN, "K", 1.);
create_temperature_unit!(DECIKELVIN, "dK", 1e-1);
create_temperature_unit!(CENTIKELVIN, "cK", 1e-2);
create_temperature_unit!(MILLIKELVIN, "mK", 1e-3);
create_temperature_unit!(MICROKELVIN, "μK", 1e-6);
create_temperature_unit!(NANOKELVIN, "nK", 1e-9);
create_temperature_unit!(PICOKELVIN, "pK", 1e-12);
create_temperature_unit!(FEMTOKELVIN, "fK", 1e-15);
create_temperature_unit!(ATTOKELVIN, "aK", 1e-18);
create_temperature_unit!(ZEPTOKELVIN, "zK", 1e-21);
create_temperature_unit!(YOCTOKELVIN, "yK", 1e-24);

#[cfg(test)]
mod tests {
    use super::*;

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
        let a = vec![METER, METER, SECOND, KILOMETER];
        let b = print_unit_from_units(a);
        println!("{b}");
        panic!()
    }

    #[test]
    fn test_derived_units() {
        let x = 5. * KILOMETER;
        let y = 10. * SECOND;
        let somethign = x * y;
        let velocity = x / y;
        let velocity_mh = velocity.to(METER / HOUR);
        assert_eq!(velocity.value, 0.5);
    }
}
