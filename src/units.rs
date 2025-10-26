use phf::phf_map;

#[derive(Debug, PartialEq)]
pub enum Dimension {
    LENGTH,
    MASS,
    TIME,
    TEMPERATURE,
    COUNT,
    DIMENSIONLESS,
}

#[derive(Debug)]
pub struct Scale {
    pub symbol: &'static str,
    pub factor: f64,
}

static SCALES: phf::Map<&'static str, Scale> = phf_map! {
    "EXA" => Scale {symbol: "E" , factor: 1e18},
    "PETA" => Scale {symbol: "P" , factor: 1e15},
    "TERA" => Scale {symbol: "T", factor: 1e12},
    "GIGA" => Scale {symbol: "G", factor: 1e9},
    "MEGA" => Scale {symbol: "M", factor: 1e6},
    "KILO" => Scale {symbol: "k", factor: 1e3},
    "HECTO" => Scale {symbol: "h", factor: 1e2},
    "DECA" => Scale {symbol: "da", factor: 1e1},
    "DECI" => Scale {symbol: "d", factor: 1e-1},
    "CENTI" => Scale {symbol: "c", factor: 1e-2},
    "MILLI" => Scale {symbol: "m", factor: 1e-3},
    "MICRO" => Scale {symbol: "u", factor: 1e-6},
    "NANO" => Scale {symbol: "n", factor: 1e-9},
    "PICO" => Scale {symbol: "p", factor: 1e-12},
    "FEMTO" => Scale {symbol: "f", factor: 1e-15},
    "ATTO" => Scale {symbol: "a", factor: 1e-18},
};

#[derive(Debug)]
pub struct SIUnit {
    pub dimension: Dimension,
    pub scale: Scale,
}

#[derive(Debug)]
pub struct Quantity {
    pub value: f64,
    pub unit: SIUnit,
}

impl Quantity {
    pub fn to(&self, other_unit: SIUnit) -> Self {
        if other_unit.dimension != self.unit.dimension {
            panic!("Dimensions are not correct!")
        }
        let ratio = other_unit.scale.factor / self.unit.scale.factor;
        Self {value: self.value * ratio, unit: other_unit}
    }
}