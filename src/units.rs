use fmtastic::Superscript;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BaseDimension {
    LENGTH,
    MASS,
    TIME,
    TEMPERATURE,
    UNITLESS,
}

#[derive(Debug, Clone, Copy, PartialEq)]
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

impl Mul for BaseUnit {
    type Output = Unit;
    fn mul(self, rhs: Self) -> Self::Output {
        self.as_unit() * rhs.as_unit()
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

impl Div<BaseUnit> for BaseUnit {
    type Output = Unit;
    fn div(self, rhs: BaseUnit) -> Self::Output {
        self.as_unit() / rhs.as_unit()
    }
}

impl Div<f64> for BaseUnit {
    type Output = Quantity;
    fn div(self, rhs: f64) -> Self::Output {
        let impl_unit = ImplBaseUnit {
            base_unit: self,
            exponent: 1,
        };
        let unit = Unit {
            base_units: vec![impl_unit],
        };
        Quantity {
            unit,
            value: 1. / rhs,
        }
    }
}
impl Div<BaseUnit> for f64 {
    type Output = Quantity;
    fn div(self, rhs: BaseUnit) -> Self::Output {
        let impl_unit = ImplBaseUnit {
            base_unit: rhs,
            exponent: -1,
        };
        let unit = Unit {
            base_units: vec![impl_unit],
        };
        Quantity { unit, value: self }
    }
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

impl Unit {
    pub fn get_units_list(&self) -> Vec<BaseUnit> {
        self.base_units.iter().map(|iu| iu.base_unit).collect()
    }
    pub fn get_unit_string(self) -> String {
        todo!()
    }
    pub fn do_dim_analysis(self) -> Vec<Dimension> {
        todo!()
    }
}

// Helper function to iterate over the implemented units and remove any that
// have zero in the exponents. If there is nothing left then it returns UNITLESS
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

impl Mul<Unit> for Unit {
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

impl Div<Unit> for Unit {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let denominator_base_units = rhs
            .base_units
            .iter()
            .map(|iu| ImplBaseUnit {
                base_unit: iu.base_unit,
                exponent: -iu.exponent,
            })
            .collect();
        self * Unit {
            base_units: denominator_base_units,
        }
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

impl UnitLike for Unit {
    fn as_unit(&self) -> Unit {
        self.clone()
    }
}

#[derive(Debug, Clone)]
pub struct Quantity {
    pub unit: Unit,
    pub value: f64,
}

impl Quantity {
    pub fn to(self, target_unit: impl UnitLike) -> Quantity {
        Quantity {
            unit: target_unit.as_unit(),
            value: self.value
                * (self.unit.calculate_conversion_factor()
                    / target_unit.calculate_conversion_factor()),
        }
    }
}

impl Add for Quantity {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let new_value = (self.value * self.unit.calculate_conversion_factor()
            + rhs.value * rhs.unit.calculate_conversion_factor())
            / self.unit.calculate_conversion_factor();
        dbg!(&new_value);
        dbg!(&rhs.value);
        dbg!(&rhs.unit.calculate_conversion_factor());
        Quantity {
            unit: self.unit.clone(),
            value: new_value,
        }
    }
}
impl Sub for Quantity {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let new_value = (self.value * self.unit.calculate_conversion_factor()
            - rhs.value * self.unit.calculate_conversion_factor())
            / self.unit.calculate_conversion_factor();
        Quantity {
            unit: self.unit.clone(),
            value: new_value,
        }
    }
}

impl Mul<Quantity> for Quantity {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl Div<Quantity> for Quantity {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        todo!()
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

create_time_unit!(MINUTE, "min", 60.);
create_time_unit!(HOUR, "hr", 3600.);

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
    fn test_adding_and_converting_units() {
        let distance_a = 5. * METER;
        let distance_b = 2. * METER;
        let meter_distance = distance_a + distance_b;
        let km_distance = meter_distance.clone().to(KILOMETER);
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
        let test_distance_meters = test_distance.clone().to(METER);
        assert_eq!(test_distance.value, 4.495);
        assert_eq!(test_distance_meters.value, 4495.);
    }

    #[test]
    fn test_printing_derived_units() {
        let unit = METER * METER * SECOND * KILOMETER;
        println!("{}", unit.get_unit_string());
        panic!()
    }

    #[test]
    fn test_derived_units() {
        let x = 5. * KILOMETER;
        let y = 10. * SECOND;
        let velocity = x / y;
        let velocity_mh = velocity.clone().to(METER / HOUR);
        assert_eq!(velocity.value, 0.5);
        assert_eq!(velocity_mh.value, 1799998.56)
    }
}
