use std::f64::consts::PI;

pub struct Parameters{
    pub x_min: f64,
    pub x_max: f64,
    pub length_x: f64,
    pub n: usize,
    pub dx: f64,
    pub t: f64, 
    pub m: usize,
    pub dt: f64,
    pub k_min: f64,
    pub k_max: f64,
    pub dk: f64,
}

impl Parameters {

    // Initiziliation of the parameters in position space.
    pub fn init(x_min: f64, x_max: f64, n: usize, t: f64, m: usize ) -> Parameters {
        let length_x: f64 = x_max - x_min;
        let dx: f64 = length_x / (n as f64);
        let dt: f64 = t / (m as f64);
        let k_min: f64 = - PI / dx;
        let k_max: f64 = - PI/ dx;
        let dk: f64 = 2.0 * PI / length_x;
        
        Parameters {
            x_min,
            x_max,
            length_x,
            n,
            dx,
            t,
            m,
            dt,
            k_min,
            k_max,
            dk
        }
    }
}

