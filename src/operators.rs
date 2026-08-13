use num_complex::Complex;

impl PsiPosition {
    pub fn hamiltonian(&mut self, potential: fn(f64) -> f64) -> Self {
        let mut x: f64 = x_min;
        for i in 0..=parameter.n {
            // multiply psi times complex exponential of position hamiltonian
            self.0[i] = self.0[i] * (Complex::new(0.0, potential(x)).exp());
            x += dx;
        }
        psi
    }
}

impl PsiMomentum {
    pub fn hamiltonian(&mut self, potential: fn(f64) -> f64) -> Self {
        let mut k: f64 = -k_range / 2.0;
        for i in 0..=parameter.n {
            // multiply psi times complex exponential of momentum hamiltonian
            psi.i = psi.i * (Complex::new(0.0, k * k).exp());
            k += dk;
        }
        psi
    }
}