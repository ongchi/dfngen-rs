use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::str::FromStr;

use itertools::zip_eq;
use parry3d_f64::na::{Point3, Vector3};
use text_io::read;
use tracing::{debug, info};

use crate::distribution::Fisher;
use crate::io::input::Input;
use crate::structures::RadiusDistribution;

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
    // {0,...}
    fn read_from_text(&mut self, file: &mut File) {
        self.extend(
            read_quoted(file, '{')
                .unwrap()
                .into_iter()
                .map(|v| v.parse().unwrap()),
        );
    }
}

impl<T> ReadFromTextFile for Vector3<T>
where
    T: FromStr + NotBool + parry3d_f64::na::Scalar + std::fmt::Debug,
    <T as FromStr>::Err: std::fmt::Debug,
{
    // {0,0,0}
    fn read_from_text(&mut self, file: &mut File) {
        let mut tmp: Vec<T> = Vec::new();
        tmp.read_from_text(file);
        if tmp.len() == 3 {
            *self = Vector3::from_row_slice(&tmp);
        } else {
            panic!("Expected 3 values in vector,")
        }
    }
}

impl ReadFromTextFile for Vec<bool> {
    // {0,...}
    fn read_from_text(&mut self, file: &mut File) {
        let mut tmp: Vec<u8> = Vec::new();
        tmp.read_from_text(file);
        self.extend(tmp.into_iter().map(|v| v != 0))
    }
}

impl ReadFromTextFile for [bool; 6] {
    // {0,0,0,0,0,0}
    fn read_from_text(&mut self, file: &mut File) {
        let mut tmp: Vec<u8> = Vec::new();
        tmp.read_from_text(file);
        for (t, s) in zip_eq(self, tmp) {
            *t = s != 0
        }
    }
}

impl<T> ReadFromTextFile for [T; 3]
where
    T: FromStr + NotBool + std::fmt::Debug,
    <T as FromStr>::Err: std::fmt::Debug,
{
    // 0
    // 0
    // 0
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
    // (0,0)
    // ...
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
    // (0,0,0)
    // ...
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

pub struct InputReader {
    file: File,
}

impl InputReader {
    pub fn new(filename: &str) -> Self {
        let file = File::open(filename).unwrap();
        Self { file }
    }
    pub fn read_value<T: ReadFromTextFile + std::fmt::Debug>(&mut self, label: &str, buf: &mut T) {
        search_var(&mut self.file, label);
        debug!("Reading from {}", label);
        buf.read_from_text(&mut self.file);
        debug!("value: {:?}", buf);
    }

    pub fn read_orien_distr(&mut self, prefix: &str) -> Vec<Fisher> {
        macro_rules! read_var {
            ($label:expr,$var_name:ident) => {
                self.read_value(&format!("{}{}:", prefix, $label), &mut $var_name);
            };
        }

        let mut angle1: Vec<f64> = Vec::new();
        let mut angle2: Vec<f64> = Vec::new();

        let mut orientation: u8 = 0;
        search_var(&mut self.file, "orientationOption:");
        orientation.read_from_text(&mut self.file);

        match orientation {
            0 => {
                read_var!("theta", angle1);
                read_var!("phi", angle2);
            }
            1 => {
                read_var!("trend", angle1);
                read_var!("plunge", angle2);
            }
            2 => {
                read_var!("dip", angle1);
                read_var!("strike", angle2);
            }
            _ => panic!("Unknown orientation option"),
        }

        let mut angle_option = false;
        search_var(&mut self.file, "angleOption:");
        angle_option.read_from_text(&mut self.file);

        if angle_option {
            let _ = angle1.iter_mut().map(|v| *v = std::f64::consts::PI / 180.);
            let _ = angle2.iter_mut().map(|v| *v = std::f64::consts::PI / 180.);
        }

        let mut kappa: Vec<f64> = Vec::new();
        read_var!("kappa", kappa);

        let mut h: f64 = 0.0;
        search_var(&mut self.file, "h:");
        h.read_from_text(&mut self.file);
        let eps = h * 1e-8;

        zip_eq(zip_eq(angle1, angle2), kappa)
            .map(|((c1, c2), k)| match orientation {
                0 => Fisher::new_with_theta_phi(c1, c2, k, eps),
                1 => Fisher::new_with_trend_plunge(c1, c2, k, eps),
                2 => Fisher::new_with_dip_strike(c1, c2, k, eps),
                _ => unreachable!(),
            })
            .collect()
    }

    pub fn read_radius_distr(&mut self, prefix: &str) -> Vec<RadiusDistribution> {
        macro_rules! read_var {
            ($label:expr,$var_name:ident) => {
                let mut $var_name: Vec<f64> = Vec::new();
                self.read_value(&format!("{}{}:", prefix, $label), &mut $var_name);
                let mut $var_name = $var_name.into_iter();
            };
        }

        macro_rules! next_val {
            ($label:expr,$var_name:ident) => {
                match $var_name.next() {
                    Some(val) => val,
                    None => panic!("Insufficient values in input file for {}{}", prefix, $label),
                }
            };
        }

        read_var!("LogMean", log_mean);
        read_var!("sd", sd);
        read_var!("LogMin", log_min);
        read_var!("LogMax", log_max);

        read_var!("alpha", alpha);
        read_var!("min", min);
        read_var!("max", max);

        read_var!("ExpMean", exp_mean);
        read_var!("ExpMin", exp_min);
        read_var!("ExpMax", exp_max);

        read_var!("const", const_);

        let mut distr: Vec<u8> = Vec::new();
        search_var(&mut self.file, &format!("{}distr:", prefix));
        distr.read_from_text(&mut self.file);

        distr
            .into_iter()
            .map(|dist_option| {
                match dist_option {
                    // Lognormal
                    1 => RadiusDistribution::new_log_normal(
                        next_val!("LogMean", log_mean),
                        next_val!("sd", sd),
                        next_val!("LogMin", log_min),
                        next_val!("LogMax", log_max),
                    ),

                    // Truncated power-law
                    2 => RadiusDistribution::new_truncated_power_law(
                        next_val!("alpha", alpha),
                        next_val!("min", min),
                        next_val!("max", max),
                    ),

                    // Exponential
                    3 => RadiusDistribution::new_exponential(
                        1. / next_val!("ExpMean", exp_mean),
                        next_val!("ExpMin", exp_min),
                        next_val!("ExpMax", exp_max),
                    ),

                    // Constant
                    4 => RadiusDistribution::new_constant(next_val!("const", const_)),

                    _ => {
                        panic!("Unknown distribution type: {}", dist_option)
                    }
                }
            })
            .collect()
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
    pub fn read_value<T: ReadFromTextFile + std::fmt::Debug>(&mut self, label: &str, buf: &mut T) {
        search_var(&mut self.file, label);
        debug!("Reading from {}", label);
        buf.read_from_text(&mut self.file);
        debug!("value: {:?}", buf);
    }

    // 0
    // ...
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

    // (0, 0)
    // ...
    pub fn read_vec_arr2(&mut self, label: &str, buf: &mut Vec<[f64; 2]>) {
        search_var(&mut self.file, label);
        debug!("Reading from {} [{}]", label, self.n_frac);
        buf.read_from_text(&mut self.file);
        debug!("value: {:?}", buf);
    }

    // (0, 0, 0)
    // ...
    pub fn read_vec_arr3(&mut self, label: &str, buf: &mut Vec<[f64; 3]>) {
        search_var(&mut self.file, label);
        debug!("Reading from {} [{}]", label, self.n_frac);
        buf.read_from_text(&mut self.file);
        debug!("value: {:?}", buf);
    }
}
