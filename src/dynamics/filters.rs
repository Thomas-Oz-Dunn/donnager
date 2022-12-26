// Filter theory

use nalgebra::*;

#[derive(Clone, Debug, PartialEq)]

pub struct KalmanFilter{
    pub syst_tm: DMatrix<f64>,
    pub obsv_tm: DMatrix<f64>,
    pub ctrl_tm: DMatrix<f64>,
    pub prcx_noise_cov: DMatrix<f64>,
    pub obsv_noise_cov: DMatrix<f64>,
    pub syst_noise_cov: DMatrix<f64>,
    pub system: DMatrix<f64>
}

impl KalmanFilter{
    
    // Create a new KF
    pub fn new(
        syst_tm: DMatrix<f64>,
        obsv_tm: DMatrix<f64>,
        ctrl_tm: DMatrix<f64>,
        prcx_noise_cov: DMatrix<f64>,
        obsv_noise_cov: DMatrix<f64>,
        syst_noise_cov: DMatrix<f64>,
        system_0: DMatrix<f64>
    ) -> KalmanFilter {
        KalmanFilter {
            syst_tm: (syst_tm), 
            obsv_tm: (obsv_tm), 
            ctrl_tm: (ctrl_tm), 
            prcx_noise_cov: (prcx_noise_cov), 
            obsv_noise_cov: (obsv_noise_cov), 
            syst_noise_cov: (syst_noise_cov), 
            system: (system_0) }
    }

    // Predict step
    pub fn predict(&mut self, control: DVector<f64>) -> DMatrix<f64> {

        // System (x) := stm * x + ctm * controls
        let system: DMatrix<f64> = 
            &self.syst_tm * self.system.transpose() + 
            &self.ctrl_tm * control.transpose();
        self.system = system.clone();

        // Posteriori noise covariance (P) := inner(stm, P) * stm + obs_noise_cov
        self.syst_noise_cov = 
            (&self.syst_tm * self.syst_noise_cov.transpose()) * &self.syst_tm 
            + &self.obsv_noise_cov;
        return system
    }

    // Update step
    pub fn update(&mut self, measurement: DVector<f64>) {

        // y = observed - inner(otm, system)
        let innovation: DVector<f64> = 
            measurement - &self.obsv_tm * self.system.transpose();
        
        // system noise = obsv_noise_cov + otm * prcx noise * otm 
        let innovation_cov: DMatrix<f64> = 
            &self.obsv_noise_cov + &self.obsv_tm * 
            (&self.prcx_noise_cov * &self.obsv_tm.transpose());

        let kalman_gain: DMatrix<f64> = 
            (&self.prcx_noise_cov * &self.obsv_tm.transpose()) * 
            innovation_cov.try_inverse().unwrap();

        // x :=  x + inner(K, y)
        let system: DMatrix<f64> = &self.system + &kalman_gain * innovation.transpose();
        self.system = system.clone();

        let dims: (usize, usize) = self.system.shape();
        let identity: DMatrix<f64> = 
            DMatrix::<f64>::identity(dims.0, dims.1);
        let z: DMatrix<f64> = 
            &identity - &kalman_gain * self.obsv_tm.transpose();

        self.syst_noise_cov = 
            (&z * self.syst_noise_cov.transpose()) * &z + 
            (&kalman_gain * self.obsv_noise_cov.transpose()) * &kalman_gain;
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::filters;

    #[test]
    fn test_kf() {

        let sys_tm: Matrix3<f64> = Matrix3::<f64>::new(
            1., 0., 0., 
            0., 1., 0., 
            0., 0., 1.);

        
        let KF: KalmanFilter = filters::KalmanFilter::new(
            sys_tm,

        );

        let result = ;
        assert_eq!(result, 4);
    }
}
