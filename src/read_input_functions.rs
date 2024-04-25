use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::str::FromStr;
use std::time::{SystemTime, UNIX_EPOCH};

use text_io::read;

use crate::read_input::Input;
use crate::structures::Point;

// // /*****************************************************************/
// // Gets multiple arrays from input/ Assumes arrays are format: {x,y,z}
// // Reads a 2D array in a 1D format.
// // Arg 1: if stream object
// // Arg 2: OUTPUT, array to place read values into
// // Arg 3: Number of rows of array we are reading
// pub fn get2dAry(file: &mut File, var: &mut [f64], row_size: usize) {
//     let mut bytes = file.bytes().map(|ch| ch.unwrap());
//
//     for _ in 0..row_size * 3 {
//         let val = read!("{}", bytes);
//         var.push(val)
//     }
// }

// // *****************************************************************
// // Gets multiple arrays from input/ Assumes arrays are format: {x,y}
// // Reads a 2D array in a 1D format.
// // Arg 1: if stream object
// // Arg 2: OUTPUT, array to place read values into
// // Arg 3: Number of rows of array we are reading
// pub fn get2dAry2(file: &mut File, var: &mut [f64], row_size: usize) {
//     let mut bytes = file.bytes().map(|ch| ch.unwrap());
//
//     for _ in 0..row_size * 2 {
//         let val = read!("{}", bytes);
//         var.push(val)
//     }
// }

// // *****************************************************************
// // Used to read in 1d arrays from input file with n Elements
// // Arg 1: ifstream object
// // Arg 2: OUTPUT, array to place read values into
// // Arg 3: Number of elements to read
// pub fn getInputAry(file: &mut File, var: &mut [f64], n_elements: usize) {
//     let mut bytes = file.bytes().map(|ch| ch.unwrap());
//
//     for _ in 0..n_elements {
//         let val = read!("{}", bytes);
//         var.push(val)
//     }
// }

// // *****************************************************************
// // Prints array of size 'nElements'
// // Arg 1: Pointer to array
// // Arg 2: Variable name
// // ARg 3: Number of elements
// pub fn printAry(var: &[f64], var_name: &str, n_elements: usize) {
//     println!("{} = {:?}", var_name, var);
// }

// // *****************************************************************
// // Prints 1-d array in 2-d array form.
// // Assumes 3 col per row
// // Arg 1: Pointer to array
// // Arg 2: Variable name
// // ARg 3: Number of row
// pub fn print2dAry(var: &[f64], var_name: &str, row_size: usize) {
//     println!("{}", var_name);
//     for i in (0..row_size * 3).step_by(3) {
//         println!("{:?}", &var[i..i + 3]);
//     }
// }

// // *****************************************************************
// // Read list of elements from file seperated by spaces
// // Arg 1: ifstream object
// // Arg 2: OUTPUT, Pointer to arary to store elements read
// // Arg 3: Number of elements to read
// pub fn getElements(file: &mut File, var: &mut [f64], n_elements: usize) {
//     let mut bytes = file.bytes().map(|ch| ch.unwrap());
//
//     for _ in 0..n_elements {
//         let val = read!("{}", bytes);
//         var.push(val);
//     }
// }

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
            println!("Variable not found: {}", search);
            panic!()
        }
    }
}

// // *******************************************************************
// // *******************************************************************
// // Checks file for being opened correectly with error msg
// // Arg 1: ifstream file object
// // Arg 2: Filename. Used for error print if there is an error
// pub fn checkIfOpen(file: &File, file_name: &str) {}

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

// // *******************************************************************
// // *******************************************************************
// // Prints rectangular coordinates. Useful for debugging.
// // Arg 1: Array that stored all rectangular coordinates.
// // Arg 2: Variable name
// // Arg 3: Number of rectangles
// fn printRectCoords(var: &[f64], varName: &str, nRectangles: usize) {
//     println!("{}:", varName);
//     for i in (0..nRectangles * 12).step_by(12) {
//         println!("{:?}", &var[i..i + 12]);
//     }
// }

// *******************************************************************
// *******************************************************************
// Gets time based seed
// Return: Seed based on the system clock
pub fn get_time_based_seed() -> u64 {
    let now = SystemTime::now();
    now.duration_since(UNIX_EPOCH).map(|t| t.as_secs()).unwrap()
}

// // checkLineLength() ********************************************************************
// // Checks if parsed line from file is the expected length.
// // Arg 1: parsed line of string (typically from splitOnWhiteSpace)
// // Arg 2: length. Expected length of line.
// // Arg 3: Filename. Used for error print if there is an error
// pub fn checkLineLength(parsed_line: &[&str], line: &str, length: usize, file_name: &str) {
//     if parsed_line.len() != length {
//         println!("Error reading {} file.", file_name);
//         println!("Line");
//         println!("{}", line);
//         println!(
//             "has a length of {}. Expected length of {}.",
//             parsed_line.len(),
//             length
//         );
//         println!("Exiting");
//         panic!()
//     }
// }

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
        global.domainVertices.push(Point {
            x: parsed_line[0].parse().unwrap(),
            y: parsed_line[1].parse().unwrap(),
            z: 0., // FIXME: uninitialized
        });
    }

    println!("Reading in Vertices Complete");
}

pub fn read_var<T>(file: &mut File) -> T
where
    T: FromStr,
    <T as FromStr>::Err: std::fmt::Debug,
{
    let mut bytes = file.bytes().map(|ch| ch.unwrap());
    let tmp: String = read!("{}", bytes);
    tmp.parse::<T>().unwrap()
}

pub fn read_vec<T>(file: &mut File) -> Vec<T>
where
    T: FromStr,
    <T as FromStr>::Err: std::fmt::Debug,
{
    let mut bytes = file.bytes().map(|ch| ch.unwrap());
    let tmp: String = read!("{}", bytes);
    tmp.strip_prefix('{')
        .map(|s| {
            s.strip_suffix('}').map(|s| {
                s.split(',')
                    .map(|s| s.parse::<T>().unwrap())
                    .collect::<Vec<T>>()
            })
        })
        .unwrap()
        .unwrap()
}

pub fn read_bool(file: &mut File) -> bool {
    read_var::<i8>(file) != 0
}

pub fn read_vec_bool(file: &mut File) -> Vec<bool> {
    let mut bytes = file.bytes().map(|ch| ch.unwrap());
    let tmp: String = read!("{}", bytes);
    tmp.strip_prefix('{')
        .map(|s| {
            s.strip_suffix('}').map(|s| {
                s.split(',')
                    .map(|s| s.parse::<i8>().unwrap() != 0)
                    .collect::<Vec<bool>>()
            })
        })
        .unwrap()
        .unwrap()
}
