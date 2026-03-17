use std::collections::HashMap;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub enum BaseDimension {
    LENGTH,
    MASS,
    TIME,
    TEMPERATURE,
    COUNT,
    DIMENSIONLESS,
}

#[derive(Debug, Clone, PartialEq)]
struct Dimension {
    dims: HashMap<BaseDimension, i32>,
}

impl Mul for Dimension {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let mut combined = self.dims.clone();
        for (k, v) in other.dims.iter() {
            let current_value = self.dims.get(k).unwrap_or(&0);
            combined.insert(*k, v + current_value);
        }
        dbg!(&combined);
        combined.retain(|_, v| *v != 0);
        Dimension { dims: combined }
    }
}
impl Div for Dimension {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let mut combined = self.dims.clone();
        for (k, v) in other.dims.iter() {
            let current_value = self.dims.get(k).unwrap_or(&0);
            combined.insert(*k, current_value - v);
        }
        dbg!(&combined);
        combined.retain(|_, v| *v != 0);
        Dimension { dims: combined }
    }
}
impl Add for Dimension {
    type Output = Result<Self, &'static str>;
    fn add(self, other: Self) -> Self::Output {
        if self.dims == other.dims {
            Ok(Dimension { dims: self.dims })
        } else {
            Err("Dimensions are not the same and can't be added.")
        }
    }
}
impl Sub for Dimension {
    type Output = Result<Self, &'static str>;
    fn sub(self, other: Self) -> Self::Output {
        if self.dims == other.dims {
            Ok(Dimension { dims: self.dims })
        } else {
            Err("Dimensions are not the same and can't be subtracted.")
        }
    }
}

/// BaseUnits are units that have a single dimension
#[derive(PartialEq, Debug, Clone)]
struct BaseUnit {
    symbol: &'static str,
    conversion_factor: f64,
    dimension: BaseDimension,
    value: f64,
}

impl Add for BaseUnit {
    type Output = Result<Self, &'static str>;
    fn add(self, other: Self) -> Self::Output {
        if self.dimension != other.dimension {
            Err("Units do not have equivalent dimentions and can't be added.")
        } else {
            let new_value = ((self.value * self.conversion_factor)
                + (other.value * other.conversion_factor))
                / self.conversion_factor;
            Ok(Self {
                symbol: self.symbol,
                conversion_factor: self.conversion_factor,
                dimension: self.dimension,
                value: new_value,
            })
        }
    }
}
impl Sub for BaseUnit {
    type Output = Result<Self, &'static str>;
    fn sub(self, other: Self) -> Self::Output {
        if self.dimension != other.dimension {
            Err("Units do not have equivalent dimentions and can't be added.")
        } else {
            let new_value = ((self.value * self.conversion_factor)
                - (other.value * other.conversion_factor))
                / self.conversion_factor;
            Ok(Self {
                symbol: self.symbol,
                conversion_factor: self.conversion_factor,
                dimension: self.dimension,
                value: new_value,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn multiplying_two_dimensions() {
        let a = Dimension {
            dims: HashMap::from([(BaseDimension::LENGTH, 2), (BaseDimension::TIME, -1)]),
        };
        let b = Dimension {
            dims: HashMap::from([(BaseDimension::LENGTH, 2), (BaseDimension::TIME, -1)]),
        };
        let c = a * b;
        assert_eq!(c.dims.len(), 2);
        assert_eq!(*c.dims.get(&BaseDimension::LENGTH).unwrap(), 4);
        assert_eq!(*c.dims.get(&BaseDimension::TIME).unwrap(), -2);
    }

    #[test]
    fn multiplying_and_cancelling() {
        // should remove factors that cancel out
        let a = Dimension {
            dims: HashMap::from([(BaseDimension::LENGTH, 2), (BaseDimension::TIME, -1)]),
        };
        let b = Dimension {
            dims: HashMap::from([(BaseDimension::LENGTH, 2), (BaseDimension::TIME, 1)]),
        };
        let c = a * b;
        assert_eq!(c.dims.len(), 1);
        assert_eq!(*c.dims.get(&BaseDimension::LENGTH).unwrap(), 4);
        assert_eq!(c.dims.get(&BaseDimension::TIME), None);
    }
    #[test]
    fn dividing_dimensions() {
        let a = Dimension {
            dims: HashMap::from([(BaseDimension::LENGTH, 2), (BaseDimension::TIME, -1)]),
        };
        let b = Dimension {
            dims: HashMap::from([(BaseDimension::LENGTH, 2)]),
        };
        let c = b / a;
        assert_eq!(c.dims.len(), 1);
        assert_eq!(*c.dims.get(&BaseDimension::TIME).unwrap(), 1);
        assert_eq!(c.dims.get(&BaseDimension::LENGTH), None);
    }
    #[test]
    fn add_and_subtract_dimensions() {
        let a = Dimension {
            dims: HashMap::from([(BaseDimension::LENGTH, 2), (BaseDimension::TIME, -1)]),
        };
        let b = Dimension {
            dims: HashMap::from([(BaseDimension::LENGTH, 2)]),
        };
        let c = a.clone() + b.clone();
        let d = a - b;
        assert!(c.is_err());
        assert!(d.is_err());
    }
    #[test]
    fn add_and_subtract_two_base_units() {
        let a = BaseUnit {
            symbol: "m",
            conversion_factor: 1.,
            dimension: BaseDimension::LENGTH,
            value: 1.,
        };
        let b = BaseUnit {
            symbol: "km",
            conversion_factor: 1000.,
            dimension: BaseDimension::LENGTH,
            value: 1.,
        };
        let c = a.clone() + b.clone();
        assert_eq!(c.clone().unwrap().value, 1001.0);
        assert_eq!(c.clone().unwrap().symbol, "m");
        assert_eq!(c.clone().unwrap().conversion_factor, 1.);
        assert_eq!(c.clone().unwrap().dimension, a.dimension);
        assert!(c.clone().is_ok());

        let c = a.clone() - b.clone();
        assert_eq!(c.clone().unwrap().value, -999.0);
        assert_eq!(c.clone().unwrap().symbol, "m");
        assert_eq!(c.clone().unwrap().conversion_factor, 1.);
        assert_eq!(c.clone().unwrap().dimension, a.dimension);
        assert!(c.clone().is_ok());
    }
    #[test]
    fn add_subtract_non_compatiable_base_units() {
        let a = BaseUnit {
            symbol: "m",
            conversion_factor: 1.,
            dimension: BaseDimension::LENGTH,
            value: 1.,
        };
        let b = BaseUnit {
            symbol: "s",
            conversion_factor: 1000.,
            dimension: BaseDimension::TIME,
            value: 1.,
        };
        let c = a.clone() - b.clone();
        let d = a.clone() + b.clone();
        assert!(c.is_err());
        assert!(d.is_err());
    }
}
