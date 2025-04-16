use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::str::FromStr;

use parry3d_f64::na::Point3;
use text_io::read;
use tracing::{debug, info};

use crate::io::input::Input;

/// Searches for variable in files, moves file pointer to position
/// after word. Used to read in varlable values
///
/// # Arguments
///
/// * `file` - File object
/// * `search` - Word to search for
pub fn search_var(file: &mut File, search: &str) {
    file.seek(SeekFrom::Start(0)).unwrap();

    let mut bytes = file.bytes().map(|ch| ch.unwrap());

    loop {
        let word: String = read!("{}", &mut bytes);
        if word == search {
            break;
        } else if word.is_empty() {
            info!(
                "Variable not found: {}",
                search.strip_suffix(":").unwrap_or(search)
            );
            std::process::exit(1)
        }
    }
}

/// Used to read in rectangualr coordinates when the user is using
/// user rectangles defined by coordinates option.
///
/// # Arguments
///
/// * `file` - file object
/// * `var` - OUTPUT. Pointer to array to store the coordinates
/// * `n_rectangles` - Number of rectangles
pub fn get_rect_coords(file: &mut File, var: &mut Vec<f64>, n_rectangles: usize) {
    let mut bytes = file.bytes().map(|ch| ch.unwrap());

    for _ in 0..n_rectangles * 12 {
        let val = read!("{}", bytes);
        var.push(val)
    }
}

/// Used to read in ellipse coordinates when the user is using
/// user ellipses defined by coordinates option.
///
/// # Arguments
///
/// * `file` - file object
/// * `out_ary` - Array to store the coordinates
/// * `n_poly` - Number of ellipses
/// * `n_vertices` - Number of points per ellipse
pub fn get_cords(file: &mut File, out_ary: &mut Vec<f64>, n_poly: usize, n_vertices: usize) {
    let mut bytes = file.bytes().map(|ch| ch.unwrap());
    let size = n_poly * n_vertices;

    for _ in 0..size * 3 {
        let val = read!("{}", bytes);
        out_ary.push(val)
    }
}

pub fn read_domain_vertices(global: &mut Input, filename: &str) {
    info!("Reading in Domain Vertices from {}", filename);

    let vertices_file = File::open(filename).unwrap();

    let reader = BufReader::new(vertices_file);
    let mut lines = reader.lines();
    let line = lines.next().unwrap().unwrap();
    let parsed_line: Vec<&str> = line.as_str().split_whitespace().collect();

    // get number of cells in x and y
    let num_of_domain_vertices: usize = parsed_line[0].parse().unwrap();
    info!(
        "There are {} Vertices on the boundary",
        num_of_domain_vertices
    );
    global.domainVertices.reserve(num_of_domain_vertices);

    for _ in 0..num_of_domain_vertices {
        let line = lines.next().unwrap().unwrap();
        let parsed_line: Vec<&str> = line.as_str().split_whitespace().collect();
        global.domainVertices.push(Point3::new(
            parsed_line[0].parse().unwrap(),
            parsed_line[1].parse().unwrap(),
            0., // FIXME: uninitialized
        ));
    }

    info!("Reading in Vertices Complete");
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

fn read_quoted(file: &mut File, quote: char) -> Result<Vec<String>, String> {
    let (left, right) = match quote {
        '(' => ('(', ')'),
        '{' => ('{', '}'),
        _ => return Err("Invalid quote character".to_string()),
    };
    let mut bytes = file.bytes().map(|c| c.unwrap());
    let mut buf = String::new();

    loop {
        let segment: String = read!("{}", bytes);
        if buf.is_empty() {
            if segment.contains(left) {
                buf.push_str(&segment);
            } else {
                return Err(format!("Expected opening quote: {}", left));
            }
            if segment.contains(right) {
                break;
            }
        } else {
            buf.push_str(&segment);
            if segment.contains(right) {
                break;
            }
        }
    }

    let lidx = buf.find(left).unwrap();
    let ridx = buf.rfind(right).unwrap();

    if lidx + 1 >= ridx {
        Ok(Vec::new())
    } else {
        Ok(buf[lidx + 1..ridx]
            .split(',')
            .map(|s| s.trim().to_string())
            .collect())
    }
}

impl<T> ReadFromTextFile for Vec<T>
where
    T: FromStr + NotBool + std::fmt::Debug,
    <T as FromStr>::Err: std::fmt::Debug,
{
    fn read_from_text(&mut self, file: &mut File) {
        self.extend(
            read_quoted(file, '{')
                .unwrap()
                .into_iter()
                .map(|v| v.parse().unwrap()),
        );
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
    T: FromStr + NotBool + std::fmt::Debug,
    <T as FromStr>::Err: std::fmt::Debug,
{
    fn read_from_text(&mut self, file: &mut File) {
        for v in self.iter_mut() {
            v.read_from_text(file);
        }
    }
}

impl<T> ReadFromTextFile for Vec<[T; 2]>
where
    T: FromStr + NotBool + std::fmt::Debug,
    <T as FromStr>::Err: std::fmt::Debug,
{
    fn read_from_text(&mut self, file: &mut File) {
        while let Ok(tmp) = read_quoted(file, '(') {
            if tmp.len() == 2 {
                self.push([tmp[0].parse().unwrap(), tmp[1].parse().unwrap()])
            } else {
                panic!("Expected 2 values in tuple,")
            }
        }
    }
}

impl<T> ReadFromTextFile for Vec<[T; 3]>
where
    T: FromStr + std::fmt::Debug,
    <T as FromStr>::Err: std::fmt::Debug,
{
    fn read_from_text(&mut self, file: &mut File) {
        while let Ok(tmp) = read_quoted(file, '(') {
            if tmp.len() == 3 {
                self.push([
                    tmp[0].parse().unwrap(),
                    tmp[1].parse().unwrap(),
                    tmp[2].parse().unwrap(),
                ])
            } else {
                panic!("Expected 3 values in tuple,")
            }
        }
    }
}

pub struct UserFractureReader {
    file: File,
    pub n_frac: usize,
}

impl UserFractureReader {
    pub fn new(filename: &str, n_frac_label: &str) -> Self {
        let mut file = File::open(filename).unwrap();

        search_var(&mut file, n_frac_label);
        let mut n_frac: usize = 0;
        n_frac.read_from_text(&mut file);

        Self { file, n_frac }
    }

    /// Read a single value from user defined fractures file
    pub fn read_value<T: ReadFromTextFile + Display>(&mut self, label: &str, buf: &mut T) {
        search_var(&mut self.file, label);
        debug!("Reading from {}", label);
        buf.read_from_text(&mut self.file);
        debug!("value: {}", buf);
    }

    pub fn read_vec<T: ReadFromTextFile + Clone + std::fmt::Debug + Display + Default>(
        &mut self,
        label: &str,
        buf: &mut Vec<T>,
    ) {
        search_var(&mut self.file, label);
        debug!("Reading from {} [{}]", label, self.n_frac);
        let mut tmp: T = Default::default();
        for _ in 0..self.n_frac {
            tmp.read_from_text(&mut self.file);
            buf.push(tmp.clone());
        }
        debug!("value: {:?}", &buf);
    }

    pub fn read_vec2(&mut self, label: &str, buf: &mut Vec<[f64; 2]>) {
        search_var(&mut self.file, label);
        debug!("Reading from {} [{}]", label, self.n_frac);
        buf.read_from_text(&mut self.file);
        debug!("value: {:?}", buf);
    }

    pub fn read_vec3(&mut self, label: &str, buf: &mut Vec<[f64; 3]>) {
        search_var(&mut self.file, label);
        debug!("Reading from {} [{}]", label, self.n_frac);
        buf.read_from_text(&mut self.file);
        debug!("value: {:?}", buf);
    }
}
