use rustfft::{FftPlanner, num_complex::Complex};
use std::f32::consts::*;


PsiPosition(Vec<Complex>)
PsiMomentum(Vec)

pub struct Parameters {
    x_range: f64,
    n: usize,
    dx: f64,
    t: f64,
    m: usize
    dt: f64,
    k: f64,
    dk: f64,
    x_min: f64
}

impl Parameters {
    pub fn init(x_range: f64, n: usize, t: f64, m: usize ) -> Parameters {
        let dx: f64 = x / (n as f64);
        let dt: f64 = t / (m as f64);
        let dk: f64 = 1.0 / x;
        let x_min: f64 = -x_range / 2.0

        Parameters {
            x_range,
            n,
            dx,
            t,
            m,
            dt,
            k_range: f64,
            dk,
            x_min
        };
    }
}

impl PsiPosition {
    pub fn init(f: fn(f64) -> f64, parameter: Parameters) -> PsiPosition {
        let mut psi: PsiPosition = // empty vec of length n;
        let mut x: f64 = x_min;

        for i in range(1, parameter.n) {
            psi.push(f(x))
            x += dx;
        }

        psi;
    }

    pub fn hamiltonian(&mut self, potential: fn(f64) -> f64) -> Self {
        let mut x: f64 = x_min;

        for i in range(1, parameter.n) {
            // multiply psi times complex exponential of position hamiltonian
            psi.i = psi.i * (Complex::new(0.0, potential(x)).exp());
            x += dx;
        }

        psi
    }

    pub fn fft(mut self) -> PsiMomentum {
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(N);
        fft.process(&self.state);
    }
}

impl PsiMomentum {
    pub fn hamiltonian(&mut self, potential: fn(f64) -> f64) -> Self {
        let mut k: f64 = -k_range / 2.0;

        for i in range(1, parameter.n) {
            // multiply psi times complex exponential of momentum hamiltonian
            psi.i = psi.i * (Complex::new(0.0, k * k).exp());
            k += dk;
        }

        psi
    }

    pub fn ifft(self) -> PsiPosition {
        let mut planner = FftPlanner::new();
        let ifft = planner.plan_fft_inverse(N);
        ifft.process(&mut self.state);
    }
}

// need a plotting library that supports complex numbers

fn main() {
    /*
    init -> loop{ plot -> position_hamiltonian -> fft -> momentum hamiltonian -> ifft}
    */
}