// Filter theory

use nalgebra::*;

#[derive(Clone, Debug, PartialEq)]

pub struct KalmanFilter{
    pub sys_trans: DMatrix<f64>,
    pub observ_trans: DMatrix<f64>,
    pub cntrl_trans: DMatrix<f64>,
    pub prcx_nosie_cov: DMatrix<f64>,
    pub obserb_noise_cov: DMatrix<f64>,
    pub system: DMatrix<f64>,
    pub post_noise_cov: DMatrix<f64>
}

impl KalmanFilter{
    
    pub fn predict(&mut self, control: DVector<f64>) -> DMatrix<f64> {
        let system: DMatrix<f64> = &self.sys_trans * self.system.transpose() + 
            &self.cntrl_trans * control.transpose();
        self.system = system.clone();
        self.post_noise_cov = 
            (&self.sys_trans * self.post_noise_cov.transpose()) * &self.sys_trans 
            + &self.obserb_noise_cov;
        return system
    }

    pub fn update(&self, measurement: DVector<f64>) {
        let y: DVector<f64> = 
            measurement - &self.observ_trans * self.system.transpose();
        let S: DMatrix<f64> = 
            &self.obserb_noise_cov + self.observ_trans * 
            (self.prcx_nosie_cov * self.observ_trans).transpose();

        let K: DMatrix<f64> = 
            (self.prcx_nosie_cov * self.observ_trans) * S.try_inverse().unwrap().transpose();
        self.system = self.system + K * y.transpose();

        let dims: (usize, usize) = self.system.shape();
        let iden: DMatrix<f64> = 
            DMatrix::<f64>::identity(dims.0, dims.1);

        self.post_noise_cov = 
            np.dot(
                np.dot(
                    iden - np.dot(
                        K, 
                        self.observ_trans), 
                    self.P), 
        	    (iden - np.dot(
                    K, 
                    self.observ_trans)).T) 
            + np.dot(
                np.dot(
                    K, 
                    self.R), 
                K.T)
    }

}