use rustfft::{FftPlanner, num_complex::Complex};
use std::f32::consts::*;

pub struct PsiPosition {
    psi: Vec<Complex<f64>>,
    dx: f64,
}

pub struct PsiMomentum {
    psi: Vec<Complex<f64>>,
    dk: f64,
}
impl PsiPosition {
    pub fn initialize(&mut self, f: fn(f64) -> f64) {
        let mut i: usize = 0;
        let mut x: f64 = self.n as f64;
        loop {
            self.psi.push(f(x)) ;
            println!("{}", f(x));
            x += L / N;
            i += 1;

            if i == N {
                break;
            }
        }
    }
    // hamiltonian operator should act on only on element of the wavefunction
    pub fn hamiltonian_position(&mut self) -> Self {
        let harmonic = (Complex::new(0.0, - M * ω * ω * (n * DX) * (n * DX) * DT)).exp();
        for n in range0..N {
            self.0[n] = self.0[n] * harmonic
        }
        
    }
    // fourier transform should map N-lattice in position space to N-lattice in momentum space

    pub fn fft(mut self) -> PsiMomentum {
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(N);
        fft.process(&self.state);
    } 

}

impl PsiMomentum {
    pub fn hamiltonian_momentum(&mut self) -> Self {
        momentum = (Complex::new(0.0, k * k / (2 * M  * H)));
        for n in range 0..P {
            self.state = self.state * momentum
            return self.state
        }
    }
    pub fn ifft(self) -> PsiPosition {
        let mut planner = FftPlanner::new();
        let ifft = planner.plan_fft_inverse(N);
        ifft.process(&mut self.state);
    }
}

fn main() {

}








