use num_complex::Complex;

use crate::parameters::*;

// PsiPosition carries the state of the wavefunciton in position space.
#[derive(Clone)]
pub struct PsiPosition(pub Vec<Complex<f64>>);

// PsiMomentum carries the state of the wavefunction in momentum space.
#[derive(Clone)]
pub struct PsiMomentum(pub Vec<Complex<f64>>);


impl PsiPosition {

    // Initialzies the wavefunction in position space given the parameters
    // and a function f which serves as the initial state of the wavefunction.
    pub fn init(f: fn(f64) -> f64, parameter: &Parameters) -> PsiPosition {
        let n: usize = parameter.n;
        let x_min: f64 = parameter.x_min;
        let dx: f64 = parameter.dx;

        let mut psi: Vec<Complex<f64>> = Vec::with_capacity(n);
        let mut x: f64 = x_min;

        for i in 0..n {
            psi.push(Complex::new(f(x), 0.0));
            x += dx;
        }
        return PsiPosition(psi);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    // Checks that the initialization function performs as expected by sampling a 
    // quadratic f(x) = x^2 over the interval[-2, 2] with N = 4. 
    // The expected values are f(-2) = 4, f(-1) = 1, f(1) = 1, f(2) = 4.
    #[test]
    fn test_psi_init() {
        fn quadratic(x: f64) -> f64 {
            x * x
        }
        let test_function: fn(f64) -> f64 = quadratic;
        let parameters: Parameters = Parameters::init(-2.0, 2.0, 4, 1.0, 1);
        let psi: PsiPosition = PsiPosition::init(test_function, &parameters);
        assert_eq!(psi.0[0], Complex::new(4.0, 0.0));
        assert_eq!(psi.0[1], Complex::new(1.0, 0.0));
        assert_eq!(psi.0[2], Complex::new(0.0, 0.0));
        assert_eq!(psi.0[3], Complex::new(1.0, 0.0));
    }
}