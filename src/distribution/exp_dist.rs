use std::cell::RefCell;
use std::rc::Rc;

use rand::distributions::Uniform;
use rand::Rng;
use rand_mt::Mt19937GenRand64;

pub struct ExpDist {
    // Maximum input to distribution
    // .999999... before being recognized
    // as 1 by the computer. (1 returns inf).
    // Initialized during Distributions constructor.
    max_input: f64,
    generator: Rc<RefCell<Mt19937GenRand64>>,
}

impl ExpDist {
    // ************** Random Uniform Random Number Generator *********************
    // Function returns a random double on [min, max]
    // Arg 1: Minimum bound
    // Arg 2: Maximum bound
    // Return: Random double on [min, max]
    pub fn new(max_decimal: f64, generator: Rc<RefCell<Mt19937GenRand64>>) -> Self {
        Self {
            max_input: max_decimal,
            generator,
        }
    }

    // ************** Random Uniform Random Number Generator *********************
    // Function returns a random double on [min, max]
    // Arg 1: Minimum bound
    // Arg 2: Maximum bound
    // Return: Random double on [min, max]
    fn unif_random(&self, min: f64, max: f64) -> f64 {
        let distribution = Uniform::new(min, max);
        self.generator.borrow_mut().sample(distribution)
    }

    // Generates a random value from the exponental distribution between the user's
    // defined minimum and maximum range.
    //
    // minVal and maxVal are the inputs needed to produce the user's minimum and maximum
    // fracture sizes. minVal and maxVal are initialized within the Distributions constructor
    // and saved in the Shape structure. Using a uniform random variable with range
    // [minVal, maxVal], the exponential distrubition will always return a value
    // within that range.
    //
    // Arg 1: Exponential Lambda (1/mean)
    // Arg 2: Minimum input (between 0 and 1)
    // Arg 3: Maximum input (between 0 and 1)
    // Return: Random number from exponential distribution described by 'lambda'
    //         and sampled with random variable between minInput and maxVal
    pub fn get_value_by_min_max_val(&self, lambda: f64, min_val: f64, max_val: f64) -> f64 {
        // Uniform distrubution on [minVal, maxVal)
        if max_val > 1. || min_val > 1. {
            // Passing 1 into exp. distribution will reuturn inf
            panic!(
        "ERROR: Passed min, or max, input value of greater than 1 to getValue() in expDist.cpp. Input must be in [0,1] interval.\n"
            )
        }

        let rand_var = self.unif_random(min_val, max_val);

        // Using inverse CDF
        // If randVar is 1, this function returns inf
        if rand_var != 1. {
            -(1. - rand_var).ln() / lambda
        } else {
            -(1. - self.max_input).ln() / lambda
        }
    }

    // *******************  Get Max Value Possible  ******************************
    // Returns maximum possible value from distribution.
    // (What the inverse CDF returns when given 0.9999... to maximum precision.
    // Inputing 1.0 into the distribution results in 'inf'.)
    //
    // Used to print warning to user if their desired maximum is larger than the
    // machine is able to produce.
    //
    // Arg 1: Lambda
    // Return: Maximum value possible due to machine precision issues
    //         before reutrning inf
    pub fn get_max_value(&self, lambda: f64) -> f64 {
        -(1. - self.max_input).ln() / lambda
    }

    // ***************************************************************************
    // *********  Compute Input to Distribution for a Certain Output  ************
    // Computes what the input value needs to be in
    // order for the distribution to return 'output'
    // Uses the distributions CDF for computation.
    //
    // Arg 1: The desired output of from the distribution (getValue())
    // Arg 2: Lambda
    // Return: The necessary input to getValue() to produce a value of 'output' */
    pub fn compute_input(output: f64, lambda: f64) -> f64 {
        1.0 - (-lambda * output).ln()
    }
}
