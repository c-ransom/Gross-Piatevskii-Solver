pub struct Parameters {
    pub x_range: f64,
    pub n: usize,
    pub dx: f64,
    pub t: f64,
    pub m: usize,
    pub dt: f64,
    pub k: f64,
    pub dk: f64,
    pub x_min: f64
}

impl Parameters {
    pub fn init(x_range: f64, n: usize, t: f64, m: usize ) -> Parameters {
        let PI: f64 = 3.14;
        let dx: f64 = x_range / ((n as f64) - 1.0);
        let dt: f64 = t / (m as f64);
        let dk: f64 = 2.0 * PI / x_range;
        let k: f64 = 0.0;
        let x_min: f64 = -x_range / 2.0;

        Parameters {
            x_range,
            n,
            dx,
            t,
            m,
            dt,
            k,
            dk,
            x_min
        }
    }
}
/* 
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parameter_init() {
        assert_eq!(3,4);
    }

}
*/