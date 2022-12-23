// Filter theory

use nalgebra::*;

pub struct KalmanFilter{
    pub F: DMatrix<f64>,
    pub H: DMatrix<f64>,
    pub B: DMatrix<f64>,
    pub P: DMatrix<f64>,
    pub Q: DMatrix<f64>,
    pub cntrl_vec:DMatrix<f64>,
    pub x: DMatrix<f64>
}

impl KalmanFilter{

    pub fn predict(&self, control: DVector<f64>) -> DMatrix<f64> {
        self.x = self.F.dot(self.x) + self.B.dot(control);
        self.P = self.F.dot(self.P).dot(self.F.transpose()) + self.Q;
        return self.x
    }

    pub fn update(&self, measurement: DVector<f64>) {
        let y: DVector<f64> = measurement - self.H.dot(self.x);
    }

}