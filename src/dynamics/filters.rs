// Filter theory

use nalgebra::*;

pub struct KalmanFilter{
    pub stm_model: DMatrix<f64>,
    pub obs_model: DMatrix<f64>,
    pub cntrl_model: DMatrix<f64>,
    pub cov_pcx_noise:DMatrix<f64>,
    pub cov_obs_noise:DMatrix<f64>,
    pub cntrl_vec:DMatrix<f64>
}

impl KalmanFilter{

    pub fn predict(&self) {
        
    }

    pub fn update(&self, measurement: DVector<f64>) {

    }

}