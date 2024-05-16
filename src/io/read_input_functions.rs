use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::str::FromStr;
use std::time::{SystemTime, UNIX_EPOCH};

use parry3d::na::Point3;
use text_io::read;

use crate::io::input::Input;

// *******************************************************************
// *******************************************************************
// Searches for variable in files, moves file pointer to position
// after word. Used to read in varlable values
// Arg 1: ifstream file object
// Arg 2: Word to search for
pub fn search_var(file: &mut File, search: &str) {
    file.seek(SeekFrom::Start(0)).unwrap();

    let mut bytes = file.bytes().map(|ch| ch.unwrap());

    loop {
        let word: String = read!("{}", &mut bytes);
        if word == search {
            break;
        } else if word.is_empty() {
            println!("Variable not found: {}", &search[..search.len() - 1]);
            std::process::exit(1)
        }
    }
}

// **********************************************************************
// **********************************************************************
// Used to read in rectangualr coordinates when the user is using
// user rectangles defined by coordinates option.
// Arg 1: ifstream file object
// Arg 2: OUTPUT. Pointer to array to store the coordinates
// Arg 3: Number of rectangles
pub fn get_rect_coords(file: &mut File, var: &mut Vec<f64>, n_rectangles: usize) {
    let mut bytes = file.bytes().map(|ch| ch.unwrap());

    for _ in 0..n_rectangles * 12 {
        let val = read!("{}", bytes);
        var.push(val)
    }
}

// **********************************************************************
// **********************************************************************
// Used to read in ellipse coordinates when the user is using
// user ellipses defined by coordinates option.
// Arg 1: ifstream file object
// Arg 2: OUTPUT. Pointer to array to store the coordinates
// Arg 3: Number of ellipses
// Arg 4: Number of points per ellipse
pub fn get_cords(file: &mut File, out_ary: &mut Vec<f64>, n_poly: usize, n_vertices: usize) {
    let mut bytes = file.bytes().map(|ch| ch.unwrap());
    let size = n_poly * n_vertices;

    for _ in 0..size * 3 {
        let val = read!("{}", bytes);
        out_ary.push(val)
    }
}

// *******************************************************************
// *******************************************************************
// Gets time based seed
// Return: Seed based on the system clock
pub fn get_time_based_seed() -> u64 {
    let now = SystemTime::now();
    now.duration_since(UNIX_EPOCH).map(|t| t.as_secs()).unwrap()
}

// Splits a line of text on white space and returns
// a of strings with the words in the line
pub fn split_on_white_space(line: &str) -> Vec<&str> {
    line.split_whitespace().collect()
}

pub fn read_domain_vertices(global: &mut Input, filename: &str) {
    println!("Reading in Domain Vertices from {}", filename);

    let vertices_file = File::open(filename).unwrap();
    // checkIfOpen(&vertices_file, filename);

    let reader = BufReader::new(vertices_file);
    let mut lines = reader.lines();
    let line = lines.next().unwrap().unwrap();
    let parsed_line = split_on_white_space(&line);

    // get number of cells in x and y
    let num_of_domain_vertices: usize = parsed_line[0].parse().unwrap();
    println!(
        "There are {} Vertices on the boundary",
        num_of_domain_vertices
    );
    global.domainVertices.reserve(num_of_domain_vertices);

    for _ in 0..num_of_domain_vertices {
        let line = lines.next().unwrap().unwrap();
        let parsed_line = split_on_white_space(&line);
        global.domainVertices.push(Point3::new(
            parsed_line[0].parse().unwrap(),
            parsed_line[1].parse().unwrap(),
            0., // FIXME: uninitialized
        ));
    }

    println!("Reading in Vertices Complete");
}

pub trait ReadFromTextFile {
    fn read_from_text(&mut self, file: &mut File);
}

trait NotBool {}
impl NotBool for u8 {}
impl NotBool for u64 {}
impl NotBool for usize {}
impl NotBool for isize {}
impl NotBool for f64 {}

impl ReadFromTextFile for String {
    fn read_from_text(&mut self, file: &mut File) {
        let mut bytes = file.bytes().map(|ch| ch.unwrap());
        *self = read!("{}", bytes);
    }
}

impl<T> ReadFromTextFile for T
where
    T: FromStr + NotBool,
    <T as FromStr>::Err: std::fmt::Debug,
{
    fn read_from_text(&mut self, file: &mut File) {
        let mut bytes = file.bytes().map(|ch| ch.unwrap());
        let tmp: String = read!("{}", bytes);
        *self = tmp.parse::<T>().unwrap()
    }
}

impl ReadFromTextFile for bool {
    fn read_from_text(&mut self, file: &mut File) {
        let mut tmp: u8 = 0;
        tmp.read_from_text(file);
        *self = tmp != 0;
    }
}

impl<T> ReadFromTextFile for Vec<T>
where
    T: FromStr + NotBool,
    <T as FromStr>::Err: std::fmt::Debug,
{
    fn read_from_text(&mut self, file: &mut File) {
        let mut bytes = file.bytes().map(|ch| ch.unwrap());
        let tmp: String = read!("{}", bytes);
        if tmp.len() > 2 {
            self.extend(
                tmp.strip_prefix('{')
                    .map(|s| {
                        s.strip_suffix('}').map(|s| {
                            s.split(',')
                                .map(|s| s.parse::<T>().unwrap())
                                .collect::<Vec<T>>()
                        })
                    })
                    .unwrap()
                    .unwrap(),
            );
        }
    }
}

impl ReadFromTextFile for Vec<bool> {
    fn read_from_text(&mut self, file: &mut File) {
        let mut tmp: Vec<u8> = Vec::new();
        tmp.read_from_text(file);
        self.extend(tmp.into_iter().map(|v| v != 0))
    }
}

impl<T> ReadFromTextFile for [T; 3]
where
    T: FromStr + NotBool,
    <T as FromStr>::Err: std::fmt::Debug,
{
    fn read_from_text(&mut self, file: &mut File) {
        let mut tmp = Vec::new();
        tmp.read_from_text(file);
        for (i, v) in tmp.into_iter().enumerate() {
            self[i] = v;
        }
    }
}
