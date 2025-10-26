use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};
use std::iter::zip;

/// The structural representation of the current IVOA UCD standard including the qualifiying
/// flags, the ucds themselves, and descriptions of those ucds.
#[derive(Debug)]
pub struct IVOAUCDWords {
    pub letters: Vec<char>,
    pub ucds: Vec<String>,
    pub descriptions: Vec<String>,
}

/// Structure to manage the list of UCD words.
///
/// Contains a list of `primary` base ucds, `secondary` base ucds. As well as a descriptions of
/// each base ucd, and the standard capitilization.
#[derive(Debug)]
pub struct UCDWords {
    pub primary: Vec<String>,
    pub secondary: Vec<String>,
    pub descriptions: HashMap<String, String>,
    pub capitilization: HashMap<String, String>,
}

impl UCDWords {
    /// Returns true if the ucd is a primary.
    pub fn is_primary(&self, ucd: String) -> bool {
        self.primary.contains(&ucd)
    }

    /// Returns true if the ucd is secondary.
    pub fn is_secondary(&self, ucd: String) -> bool {
        self.secondary.contains(&ucd)
    }

    /// Returns the description of the given (case insensitive) ucd string.
    pub fn get_description(&self, ucd: String) -> Option<&String> {
        self.descriptions.get(&ucd)
    }

    /// Returns the standard capitilization of the given (case insensitive) ucd string.
    pub fn normalize_capitilization(&self, ucd: String) -> Option<&String> {
        self.capitilization.get(&ucd)
    }

    /// Returns true is the ucd string is a valid ucd.
    pub fn check_ucd(&self, ucd: &str) -> bool {
        let ucd_lower = ucd.to_lowercase();
        let parts: Vec<&str> = ucd_lower.split(';').collect();

        let valids: Vec<bool> = parts
            .iter()
            .enumerate()
            .map(|(i, &part)| {
                if i == 0 {
                    self.is_primary(String::from(part))
                } else {
                    self.is_secondary(String::from(part))
                }
            })
            .collect();

        !valids.contains(&false)
    }
}

/// Loads the list of ucd words.
pub fn load_ivoa_ucd_words() -> io::Result<IVOAUCDWords> {
    let mut file = File::open("data/ucd_1p6.dat")?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    let lines: Vec<&str> = contents
        .lines()
        .filter(|line| !line.starts_with("#"))
        .collect();
    let mut letters: Vec<char> = Vec::new();
    let mut ucds: Vec<String> = Vec::new();
    let mut descriptions: Vec<String> = Vec::new();
    for line in lines {
        let row: Vec<&str> = line.split('|').map(|col| col.trim()).collect();
        letters.push(row[0].chars().next().unwrap());
        ucds.push(String::from(row[1]));
        descriptions.push(String::from(row[2]))
    }
    Ok(IVOAUCDWords {
        letters,
        ucds,
        descriptions,
    })
}

/// Reads the read in UCD words that creates the UCDWords struct.
///
/// #Example
/// ```
/// use crate::astroxide::ucd_checker::{load_ivoa_ucd_words, read_ucd_words};
///
/// let words = load_ivoa_ucd_words().expect("File not found");
/// let ucd_words = read_ucd_words(words);
///
/// let valid_string = "meta.id;meta.main";
/// let also_valid_string = "pOs.EQ.Ra;mETa.maIN"; // case insensitive.
/// let not_valid = "meta.main;arith"; // arith is a primary and meta.main is secondary.
///
/// assert!(ucd_words.check_ucd(valid_string));
/// assert!(ucd_words.check_ucd(also_valid_string));
/// assert!(!ucd_words.check_ucd(not_valid));
/// ```
pub fn read_ucd_words(words: IVOAUCDWords) -> UCDWords {
    let mut primary: Vec<String> = Vec::new();
    let mut secondary: Vec<String> = Vec::new();
    let mut descriptions: HashMap<String, String> = HashMap::new();
    let mut capitilization: HashMap<String, String> = HashMap::new();
    let letters = words.letters;
    let ucds = words.ucds;
    let description = words.descriptions;
    for (letter, ucd) in zip(letters, &ucds) {
        if "QPEVC".contains(letter) {
            primary.push(ucd.clone().to_lowercase())
        }
        if "QSEVC".contains(letter) {
            secondary.push(ucd.clone().to_lowercase())
        }
    }

    for (ucd, descr) in zip(&ucds, description) {
        descriptions.insert(ucd.clone().to_lowercase(), descr);
    }

    for ucd in ucds {
        capitilization.insert(ucd.clone().to_lowercase(), ucd);
    }
    UCDWords {
        primary,
        secondary,
        descriptions,
        capitilization,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reading_in_file() {
        let ucds = load_ivoa_ucd_words().unwrap();
        assert_eq!(ucds.letters[0], 'Q');
        assert_eq!(ucds.descriptions[100], "Various instrumental parameters");
        assert_eq!(ucds.ucds.last().unwrap(), "time.start");
    }

    #[test]
    fn building_ucd_words() {
        let words = load_ivoa_ucd_words().unwrap();
        let ucd_words = read_ucd_words(words);
        assert_eq!(ucd_words.secondary.first().unwrap(), "arith");
        assert_eq!(ucd_words.secondary.get(1).unwrap(), "arith.diff");
        assert_eq!(ucd_words.primary.first().unwrap(), "arith");
        assert_eq!(ucd_words.primary.get(1).unwrap(), "arith.factor");

        assert_eq!(
            ucd_words.descriptions.get("em.ir.15-30um").unwrap(),
            "Infrared between 15 and 30 micron"
        );
        assert_eq!(
            ucd_words.capitilization.get("em.ir.15-30um").unwrap(),
            "em.IR.15-30um"
        );
    }

    #[test]
    fn test_ucd_methods() {
        let words = load_ivoa_ucd_words().unwrap();
        let ucd_words = read_ucd_words(words);
        assert!(ucd_words.is_primary(String::from("arith")));
        assert!(ucd_words.is_secondary(String::from("arith")));
        assert!(ucd_words.is_secondary(String::from("arith.diff")));
        assert!(ucd_words.is_primary(String::from("arith.factor")));
        assert!(!ucd_words.is_primary(String::from("arith.diff")));
        assert!(!ucd_words.is_secondary(String::from("arith.factor")));

        assert_eq!(
            ucd_words
                .get_description(String::from("em.ir.15-30um"))
                .unwrap(),
            "Infrared between 15 and 30 micron"
        );
        assert_eq!(
            ucd_words
                .normalize_capitilization(String::from("em.ir.15-30um"))
                .unwrap(),
            "em.IR.15-30um"
        );
    }

    #[test]
    fn test_check_ucd() {
        let valid_ucd = "arith;em.IR.8-15um";
        let not_valid1 = "em.IR.8-15um";
        let valid_case_insensitive = "arItH;EM.IR.8-15um";
        let not_valid_order = "em.IR.8-15um;arith";
        let not_valid_random = "arith;foobar";
        let not_valid_random_2 = "foobar;em.IR.8-15um";

        let words = load_ivoa_ucd_words().unwrap();
        let ucds = read_ucd_words(words);
        assert!(ucds.check_ucd(valid_ucd));
        assert!(!ucds.check_ucd(not_valid1));
        assert!(ucds.check_ucd(valid_case_insensitive));
        assert!(!ucds.check_ucd(not_valid_order));
        assert!(!ucds.check_ucd(not_valid_random));
        assert!(!ucds.check_ucd(not_valid_random_2));
    }
}
