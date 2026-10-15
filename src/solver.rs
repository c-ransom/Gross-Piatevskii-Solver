use num_complex::Complex;
use crate::wavefunction::*;

use crate::parameters::Parameters;

impl PsiPosition {

    // Acts the position piece of hamiltonian on wavefunction
    pub fn hamiltonian_position(mut self, parameter: &Parameters ) -> Self {
        let dt: f64 = parameter.dt;
        let dx: f64 = parameter.dx;
        let mut x: f64 = parameter.x_min;

        for amplituide in self.0.iter_mut() {

            // The position hamiltonian (Ĥ_x) for GPE is V(x,t) + |Ψ|^2 where V(x,t) is harmonic oscillator (quadratic).
            let hamiltonian: f64 = x * x / 2.0 + amplituide.norm();

            // The time evolution is declared by e^(iĤ_x * dt)
            let time_evolution: Complex<f64> = Complex::new(0.0, - hamiltonian * dt / 2.0).exp();
            
            // Since position hamiltonian is diagonilzed in position space
            // the hamiltonian just amounts to multplication. by e^(i )
            *amplituide = (*amplituide) * time_evolution ;
            
            x += dx;
        }
        return self;
    }

}

impl PsiMomentum {
    // Acts the kinetic piece of the hamiltonian on wavefunction
    pub fn hamiltonian_momentum(mut self, parameter: &Parameters) -> Self {
        let dt: f64 = parameter.dt;
        let dk: f64 = parameter.dk;
        let mut k: f64 = parameter.k_min;
        
        for amplituide in self.0.iter_mut() {
            let hamiltonian: f64 = k * k;
            let time_evolution: Complex<f64> = Complex::new(0.0, - hamiltonian * dt).exp();

            *amplituide = (*amplituide) * time_evolution;
            k += dk;
        }
        return self;
    }
}