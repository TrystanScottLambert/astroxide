use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};
use std::iter::zip;

#[derive(Debug)]
pub struct IVOAUCDWords {
    pub letters: Vec<char>,
    pub ucds: Vec<String>,
    pub descriptions: Vec<String>,
}

#[derive(Debug)]
pub struct UCDWords {
    pub primary: Vec<String>,
    pub secondary: Vec<String>,
    pub descriptions: HashMap<String, String>,
    pub capitilization: HashMap<String, String>,
}

impl UCDWords {
    pub fn is_primary(&self, ucd: String) -> bool {
        self.primary.contains(&ucd)
    }

    pub fn is_secondary(&self, ucd: String) -> bool {
        self.secondary.contains(&ucd)
    }

    pub fn get_description(&self, ucd: String) -> Option<&String> {
        self.descriptions.get(&ucd)
    }

    pub fn normalize_capitilization(&self, ucd: String) -> Option<&String> {
        self.capitilization.get(&ucd)
    }
}

fn load_ivoa_ucd_words() -> io::Result<IVOAUCDWords> {
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
            println!("{}", ucd.clone());
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
}
