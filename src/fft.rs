use num_complex::Complex;
use std::f64::consts::PI;
use crate::parameters::Parameters;
use crate::wavefunction::*;


impl PsiPosition {
    // This is the ordinary fourier transform algorithm O(N^2).
    pub fn ft(self, parameter: &Parameters) -> PsiMomentum {
        let mut ft_psi: Vec<Complex<f64>> = Vec::new();

        for k in 0..parameter.n {
            let mut x_k: Complex<f64> = Complex::new(0.0, 0.0);

            for n in 0..parameter.n {
                let phase: f64 = - 2.0 * PI * (k as f64) * (n as f64) / (parameter.n as f64);
                x_k += (self.0[n])
                    * (Complex::new(0.0, phase).exp()); 
            }
            ft_psi.push(x_k);
        }
        PsiMomentum(ft_psi)
    }

    // This is the fast fourier transform algorithm O(Nlog(N)).
    pub fn fft(self) -> PsiMomentum {
        let n = self.0.len();
        if n == 1 {
            // Solves the base case.
            return PsiMomentum(vec![self.0[0]]);
        } else {

        // Solves the inductive case.
            let mut even = PsiPosition(Vec::new());
            let mut odd = PsiPosition(Vec::new());

            // Splits the array into two arrays where one has the even indexed components
            // and the other has the odd indexed components each with length N/2.
            for i in 0..n {
                if i % 2 == 0 {
                    even.0.push(self.0[i]);
                } else {
                    odd.0.push(self.0[i]);
                }
            }

            let even_fft = even.fft();
            let odd_fft = odd.fft();
            
            // Calculates the twiddle factor, w = e^-(2πi/N).
            let w = Complex::new(0.0, - 2.0 * PI / n as f64).exp(); 
            let mut x = PsiMomentum(vec![
                Complex::new(0.0, 0.0);
            n]);

            
            let normalization: f64 = 2.0_f64.sqrt();

            for k in 0..(n / 2) {
                let w_k = w.powu(k as u32);
                
                // Computes the fourier transform in the recursive form.
                // Normalization ensures the wavefunction does 
                // not blow up after repeated iterations of fourier transform.
                x.0[k] = (even_fft.0[k] + w_k * odd_fft.0[k]) / normalization;
                x.0[k + n / 2] = (even_fft.0[k] - w_k * odd_fft.0[k]) / normalization;
            }
            return x
        }
    }
} 

impl PsiMomentum {
    // This is the ordinary inverse fourier transform algorithm O(N^2).
    pub fn ift(self, parameter: &Parameters) -> PsiPosition {
        let mut ift_psi = Vec::new();
        
        for n in 0..parameter.n {
            let mut x_n: Complex<f64> = Complex::new(0.0, 0.0);

            for k in 0..parameter.n {
                let phase: f64 = - 2.0 * PI * (n as f64) * (k as f64) / (parameter.n as f64);
                x_n += self.0[k] * Complex::new(0.0, phase).exp()  
            }
            ift_psi.push(x_n);
        }
        PsiPosition(ift_psi)
    }

    // This is the fast fourier transform algorithm O(Nlog(N)).
    pub fn ifft(self) -> PsiPosition {
        let n = self.0.len();
        if n == 1 {
            // Calculates the base case of the fft.
            return PsiPosition(vec![self.0[0]]);

        } else {
            
        // Solves the inductive case.
            let mut even = PsiMomentum(Vec::new());
            let mut odd = PsiMomentum(Vec::new());

            // Splits the array into two arrays where one has the even indexed components
            // and the other has the odd indexed components each with length N/2.
            for i in 0..n {
                if i % 2 == 0 {
                    even.0.push(self.0[i]);
                } else {
                    odd.0.push(self.0[i]);
                }
            }
            
            let even_ifft = even.ifft();
            let odd_ifft = odd.ifft();

            // Calculates the twiddle factor, w = e^(2πi/N).
            let w = Complex::new(0.0, 2.0 * PI / n as f64).exp(); 
            let mut x = PsiPosition(vec![
                Complex::new(0.0, 0.0);
             n]);

            let normalization: f64 = 2.0_f64.sqrt();
            for k in 0..(n / 2) {
                let w_k = w.powu(k as u32);

                // Computes the fourier transform in the recursive form.
                // Normalization ensures the wavefunction does 
                // not blow up after repeated iterations of fourier transform.
                x.0[k] = (even_ifft.0[k] + w_k * odd_ifft.0[k]) /normalization;
                x.0[k + n / 2] = (even_ifft.0[k] - w_k * odd_ifft.0[k]) / normalization;
            }
            return x
        }
    }

}
    
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    // This test checks if the fourier transform computes simple case of sin(x) from [0,2π] with N = 4
    // where we should expect peaks at frequency peaks at k = -1, 1.
    fn test_fourier() {
        let test_function: fn(f64) -> f64 = f64::sin;
        let parameters: Parameters = Parameters::init(0.0, 2.0 * PI , 4, 1.0, 1);
        let psi: PsiPosition = PsiPosition::init(test_function, &parameters);
        println!("These are the components of sin in position space: {}, {}, {}, {}", 
        psi.0[0], psi.0[1], psi.0[2], psi.0[3]);

        let ft_psi = psi.ft(&parameters);
        println!("These are the componets of sin in momentum space: {}, {}, {}, {}", 
        ft_psi.0[0],ft_psi.0[1], ft_psi.0[2], ft_psi.0[3]);

        assert!((ft_psi.0[0] - Complex::new(0.0, 0.0)).norm() < 0.001);
        assert!((ft_psi.0[1] - Complex::new(0.0,-2.0)).norm() < 0.001);
        assert!((ft_psi.0[2] - Complex::new(0.0, 0.0)).norm() < 0.001);
        assert!((ft_psi.0[3] - Complex::new(0.0, 2.0)).norm() < 0.001);

    }

    #[test]
    // This test checks if the fast fourier transform computes simple case of sin(x) from [0,2π] with N = 4
    // where we should expect peaks at frequency peaks at k = -1, 1.
    fn test_fft() {
        let test_function: fn(f64) -> f64 = f64::sin;
        let parameters: Parameters = Parameters::init(0.0, 2.0 * PI, 4, 2.0, 1);
        let psi: PsiPosition = PsiPosition::init(test_function, &parameters);
        println!("These are the components of sin in position space: {}, {}, {}, {}",
    &psi.0[0], &psi.0[1], &psi.0[2], &psi.0[3]);

    let fft_psi = psi.fft();
    println!(
        "These are the components of sin in momentum: {}, {}, {}, {}",
        fft_psi.0[0], fft_psi.0[1], fft_psi.0[2], fft_psi.0[3]);

    assert!((fft_psi.0[0] - Complex::new(0.0, 0.0)).norm() < 0.001);
    assert!((fft_psi.0[1] - Complex::new(0.0,-1.0)).norm() < 0.001);
    assert!((fft_psi.0[2] - Complex::new(0.0, 0.0)).norm() < 0.001);
    assert!((fft_psi.0[3] - Complex::new(0.0, 1.0)).norm() < 0.001);
    }

    #[test]
    // This test checks if the fft and ifft form identity for simples case of sin(x) from [0,2π] with N = 4.
    fn test_ifft_fft() {
        let test_function: fn(f64) -> f64 = f64::sin;
        let parameters: Parameters = Parameters::init(0.0, 2.0 * PI, 4, 2.0, 1);
        
        let psi: PsiPosition = PsiPosition::init(test_function, &parameters);
        let psi_clone = psi.clone();
        let ifft_fft_psi: PsiPosition = psi.fft().ifft();

    for i in 0..ifft_fft_psi.0.len() {
        let acc_error: f64 = 0.1;
        
        let exp_error: f64 = (psi_clone.0[i] - ifft_fft_psi.0[i]).norm();
        println!(
            "Component {i} error: {exp_error} \n 
            Initialied wavefunction values: {} \n
            ifft_fft wavefunction values: {} )", psi_clone.0[i], ifft_fft_psi.0[i]);

        assert!( exp_error < acc_error);
    }

    }
}
