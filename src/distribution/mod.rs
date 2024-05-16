use crate::fracture::insert_shape::{get_family_number, shape_type};
use crate::io::input::Input;
use crate::structures::Shape;
use exp_dist::ExpDist;
use rand_mt::Mt19937GenRand64;
use std::cell::RefCell;
use std::io::Read;
use std::rc::Rc;

pub mod exp_dist;
pub mod generating_points;

// The Distributions class is used to hold all the custom distribution classes
// that are used by DFNgen. As of now, the only distribution we have completely customized
// is the exponential distribution class ExpDist. The customizations allow us
// to limit the range of numbers the distribution is able to produce. This allows
// the user to define the minimum and maximum size fractures they want from the
// distribution.
pub struct Distribution {
    pub max_input: f64,
    pub exp_dist: ExpDist,
}

// The Distributions class was created specifically for
// easy use of a customized exponential distribution
// which allows users to specify the min and max limit
// which the distribution function can return.
// (exponentail distribution function truncated on both sides)
//
// The Distributions class was created keeping in mind
// other distributions may be needed later.
//
// Currently, exponential distribution is the only
// distribution contained in the Distributions class.
impl Distribution {
    // Initialize maxInput. maxInput is the maximum double
    // less than 1. (0.9999... before being recognized as 1)
    //
    // Initialize exponential distribution class.
    //
    // Arg 1: Random number generator */
    pub fn new(
        input: &Input,
        generator: Rc<RefCell<Mt19937GenRand64>>,
        shape_families: &mut [Shape],
    ) -> Self {
        // Maximum double less than 1.0
        // (0.9999... before being recognized as 1)
        let max_input = Self::get_max_decimal_for_double();
        // Init exponential dist class
        let exp_dist = ExpDist::new(max_input, generator);

        let distrubutions = Self {
            max_input,
            exp_dist,
        };

        // User input check and exp dist init
        distrubutions.check_distribution_user_input(input, shape_families);

        distrubutions
    }

    // ***********************************************************
    // ******************  Get Max Digits  ***********************
    // Calculates and returns the maximum precision, in number of
    // digits, after the decimal place for a number between 0 and 1.
    fn get_max_digits() -> u32 {
        (f64::MANTISSA_DIGITS as f64 * std::f64::consts::LOG10_2).floor() as u32 + 2
    }

    // /***********************************************************/
    // /************** Get Max Number Less Than 1.0  **************/
    // /*! Returns the largest number less than 1.0 before being
    //     considered 1.0 by the computer (0.99999...). */
    fn get_max_decimal_for_double() -> f64 {
        let mut temp = String::with_capacity(32);
        let length = Self::get_max_digits();

        temp.push('.');
        for _ in 1..length {
            temp.push('9')
        }

        temp.parse().unwrap()
    }

    // ********************************************************************************/
    // *** Initialize and  Error Check On User Input for Exponendtail Distribution  ***/
    // MANDATORY FUNCTION FOR USING EXPONENTIAL DISTRIBUTION
    // This function currently only error checks the exponential
    // distribution option. The code has been outlined to
    // add other distributions easily in the future.
    // A custom exponential distribution function was needed to allow
    // the user to define the minimum and maximum value the distribution
    // would return, without sampling randomly and throwing away numbers
    // larger than the maximum or less than the minimum user defined limit.
    //
    // This function checks that the user's min and max limit for exponential
    // distribution is acceptible. That is, that the machine is able to produce
    // the min and maximum given. For example, if the mean is set very small,
    // and the user defined maximim is set very large, it may be the case that the
    // largest number able to be produced from the distribution will be less than
    // the user defined maximum. This check will warn the user
    // if this is the case.
    //
    // This function is MANDATORY for the min and max input options to work.
    // The function initializes the minimum and maximum input to the
    // distrubution function in order to return a value which is between
    // the user's min and max distribution limit from the input file.
    //
    // Arg 1: Array of all stochastic fracture families */
    fn check_distribution_user_input(&self, global: &Input, shape_families: &mut [Shape]) {
        let mut input;
        for (i, shape) in shape_families.iter_mut().enumerate() {
            match shape.distribution_type {
                // Lognormal
                1 => {}
                // Truncated power-law
                2 => {}
                // Exponential
                3 => {
                    // Check exponential minimum value
                    input = ExpDist::compute_input(shape.exp_min, shape.exp_lambda);

                    if input >= 1. {
                        println!(
                "\n\nWARNING: The defined minimum value, {}, for {} family {} will not be able to be produced by the exponential distribution due to precision issues.", shape.exp_min, shape_type(shape), get_family_number(global, i as isize, shape.shape_family)
                    );
                        println!(
                "The minimum value is too large. The largest value the distribution can produce with current parameters is {}",
                self.exp_dist.get_max_value(shape.exp_lambda)
                    );
                        println!("Please adjust the minimum value, or the mean, and try again.");
                        panic!()
                    } else {
                        shape.min_dist_input = input;
                    }

                    // Check exponential maximum value
                    input = ExpDist::compute_input(shape.exp_max, shape.exp_lambda);

                    if input >= 1. {
                        println!(
                    "\n\nWARNING: The defined maximum value, {}, for {} family {} will not be able to be produced by the exponential distribution due to precision issues.",
                    shape.exp_max,
                    shape_type(shape),
                    get_family_number(global, i as isize, shape.shape_family)
                );
                        println!(
                    "The largest value the distribution can produce with current parameters is {}",
                self.exp_dist.get_max_value(shape.exp_lambda)
                    );
                        println!(
                "Press Enter to automatically adjust this exponential maximum value to {} (q to Quit)",
                self.exp_dist.get_max_value(shape.exp_lambda)
                    );

                        // Prompt user to press enter to continue or q to quit
                        Self::quit_or_continue();
                        let max = self.exp_dist.get_max_value(shape.exp_lambda);

                        //check that the max is not less or equal to the min
                        if max <= shape.exp_min {
                            println!(
                    "ERROR: The maximum exponetnial distribution radius possible for {} family {} is less than or equal to the minimum exponential distribution radius.",
                    shape_type(shape),
                    get_family_number(global, i as isize, shape.shape_family)
                        );
                            println!(
                    "Please adjust the exponential distribution parameters in {} family {}",
                    shape_type(shape),
                    get_family_number(global, i as isize, shape.shape_family)
                        );
                            panic!()
                        }
                    }

                    shape.max_dist_input = input;
                }
                _ => {
                    unreachable!()
                }
            }
        }
    }

    // ******************************************************
    // *************** Quit or Continue *********************
    // Stops program execution and waits for input from user
    // to continue.
    // 'q' or 'Q' to quit program and exit
    // enter/return key to continue with program execution.*/
    fn quit_or_continue() {
        loop {
            let mut input = [0u8];
            let _ = std::io::stdin().read(&mut input);
            match input[0] as char {
                'q' | 'Q' => panic!(),
                '\r' | '\n' => break,
                _ => println!("Invalid input."),
            }
        }
    }
}
