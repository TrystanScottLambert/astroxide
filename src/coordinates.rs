pub enum Dimension {
    length,
    mass,
    time,
    temperature,
    count,
    dimensionless,
}

struct Multiple {
    prefix: &str,
    symbol: &str,
    factor: f64,
}


pub enum Scale {
    exa,
    peta,
    tera,
    giga,
    mega,
    kilo,
    hecto,
    deca,
    deci,
    centi,
    milli,
    micro,
    nano,
    pico,
    femto,
    atto,
}

pub struct Unit {
    dimension: Dimension,
    scale: Scale,
}

pub struct Quantity {
    value: f64,
    unit: Unit
}