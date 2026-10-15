mod wavefunction;
mod parameters;
mod fft;
mod solver;
use std::io;
use crate::wavefunction::*;
use crate::parameters::Parameters;
use crate::fft::*;
use crate::solver::*;
use csv::*;
use std::process;
use std::f64::consts::PI;
use num_complex::Complex;
use std::error::Error;
use std::result::Result;


// Advances the wavefunction one time step fowards
// using split-operator method.

fn evolve(psi: PsiPosition, parameters: &Parameters) -> PsiPosition  {

    // Apply position operators in position space
    let psi = psi.hamiltonian_position(parameters);
    // Fourier transform into momentum space
    let psi = psi.fft();
    // Apply momentum operators in momentum space
    let psi = psi.hamiltonian_momentum(parameters);
    // Fourier transform back into position space
    let psi = psi.ifft();
    psi
}

fn main() -> Result<(), Box<dyn Error>> {

    println!("Beginning simulation");
    println!("Initializing Parameters");

    let parameters: Parameters = 
    Parameters::init(-1.0, 1.0, 2048, 10.0, 16384);
    let n: usize = parameters.n;
    let m: usize = parameters.m;

        fn gaussian(x: f64) -> f64 {
        (- (x * x)).exp()
    }

    println!("Initializing Wavefunction");
    
    // Initializes the wavefunction as a gaussian function.
    let mut psi: PsiPosition = 
        PsiPosition::init(gaussian, &parameters);
    
    // Create the CSV writer once.
    let mut writer = csv::Writer::from_path("gpe_data.csv")?;
    
    // Evolve the system through time.
    for step in 0..m {

        psi = evolve(psi, &parameters);

        // Calculate the density at each spatial grid point.
        let density: Vec<String> = psi.0
            .iter()
            .map(|z| z.norm_sqr().to_string())
            .collect();

        // Write one row for this time step.
        writer.write_record(&density)?;

        println!("Completed step {step}");
    }
    
        // Ensure buffered data is written.
        writer.flush()?;
    
        Ok(())
}


#[cfg(test)]
mod tests {

    use super::*; 

    // Checks whether the kinetic part of hamiltonian behaves as expected for simple case of sin(x) [0, 2π] with N = 4. 
    #[test]
    fn test_main() {
        
        let parameters: Parameters = 
        Parameters::init( 0.0, 2.0 * PI, 4, 4.0, 4);
        let dt = parameters.dt;

        let test_function: fn(f64) -> f64 = f64::sin;

        let mut psi: PsiPosition = 
        PsiPosition::init(test_function, &parameters);
        
        let truth: PsiPosition = 
        PsiPosition(vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, (-PI * PI / 8.0 - 1.0) / 2.0 * dt).exp(),
            Complex::new(0.0, 0.0),
          - Complex::new(0.0, (-9.0 * PI * PI / 8.0 - 1.0) / 2.0 * dt).exp()
            ]);
        
        psi = psi.hamiltonian_position(&parameters);


        for i in 0..psi.0.len() {
            let acc_error: f64 = 1.0;
            let exp_error: f64 = (psi.0[i] - truth.0[i]).norm();
            println!(
                    "Component {} error: {} \n 
                    Output value: {} \n 
                    True Value: {})",
                     i, exp_error, psi.0[i], truth.0[i]);
            assert!( exp_error < acc_error);
        }
    
    }
}
 