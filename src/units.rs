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

#[derive(PartialEq, Debug)]
struct Unit {
    symbol: &'static str,
    conversion_factor: f64,
    dimensions: Dimension,
}
impl Add for Unit {
    type Output = Result<Self, &'static str>;
    fn add(self, other: Self) -> Self::Output {
        let new_dimension = self.dimensions + other.dimensions;
        if new_dimension.is_err() {
            Err("Units do not have equivalent dimentions and can't be added.")
        } else {
            let new_string = format!("{} {}", self.symbol, other.symbol).as_str();
        }
    }
}

#[derive(PartialEq, Debug)]
pub struct Length {
    value: f64,
    unit: Unit,
}

impl Length {
    pub fn new(value: f64, unit: Unit) -> Self {
        Self { value, unit }
    }

    pub fn convert(&self, other: Unit) -> Self {
        todo!()
    }
}
impl Add for Length {
    type Output = Result<Self, &'static str>;
    fn add(self, other: Self) -> Self::Output {
        let new_unit = self.unit + other.unit;
        let new_value = self.value + other.value;
        Ok(Self {
            value: new_value,
            unit: new_unit?,
        })
    }
}

impl Sub for Length {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            value: self.value - other.value,
            symbol: self.symbol,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // #[test]
    // fn adding_two_lengths() {
    //     let a = Length::new(2.0, "km") + Length::new(2.0, "km");
    //     assert_eq!(a, Length::new(4.0, "km"));
    // }

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
}
