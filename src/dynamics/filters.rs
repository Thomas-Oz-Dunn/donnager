// Filter theory

use nalgebra::*;

pub struct KalmanFilter{
    pub stm_model: Matrix<f64>,
    pub obs_model: Matrix<f64>,
    pub cntrl_model: Matrix<f64>,
    pub cov_pcx_noise:,
    pub cov_obs_noise:,
    pub cntrl_vec:
}

impl KalmanFilter{

    pub fn predict(&self) {

    }

    pub fn update(&self, obs) {

    }

}