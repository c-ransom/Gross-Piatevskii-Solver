use num_complex::Complex;
use std::f64::consts::PI;

impl PsiPosition {
    pub fn ft(self) -> PsiMomentum {
        let mut psi = Vec::new();
        for k in 0..parameter.n {
            let mut x_k: Complex<f64> = Complex::new(0.0, 0.0);
            for n in 0..parameter.n {
                x_k += (psi[n])*(Complex::new(0.0, 2* PI * k * n / parameter.n).exp());
            }
            psi.push(x_k);
        }
    }
}

impl PsiMomentum {

}

#[cfg(test)]
mod tests {
    super::*;

    #[test]
    fn test_fourier() {
        let test_function: fn(f64) -> f64 = sin();
        let parameters: Parameters = Parameters::init(2 * PI, 4, 0, 0);
        let psi: PsiPosition = PsiPosition::init(test_function, parameters);
        let ft_psi = ft(psi);
        assert_eq!()
    }
}
