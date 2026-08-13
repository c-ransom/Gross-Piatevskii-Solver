use num_complex::Complex;
use crate::parameters::Parameters;
use std::f64::consts::PI;

struct PsiPosition(Vec<Complex<f64>>);
struct PsiMomentum(Vec<Complex<f64>>);

impl PsiPosition {
    pub fn init(f: fn(f64) -> f64, parameter: Parameters) -> PsiPosition {
        let mut psi = Vec::new();
        let mut x: f64 = parameter.x_min;

        for i in 0..parameter.n {
            psi.push(Complex::new(f(x), 0.0));
            println!("{}", psi[i]);
            x += parameter.dx;
            
        }

        PsiPosition(psi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_psi_init() {
        fn quadratic(x: f64) -> f64 {
            x * x
        }
        let test_function: fn(f64) -> f64 = quadratic;
        let parameters: Parameters = Parameters::init(2.0, 3, 1.0, 1);
        let psi: PsiPosition = PsiPosition::init(test_function, parameters);
        assert_eq!(psi.0[0], Complex::new(1.0, 0.0));
        assert_eq!(psi.0[1], Complex::new(0.0, 0.0));
        assert_eq!(psi.0[2], Complex::new(1.0, 0.0));

    }
}