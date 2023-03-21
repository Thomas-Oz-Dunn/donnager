/* 
Propogator 
*/

use nalgebra::Vector3;
use chrono::{DateTime as DateTime, Utc};

use crate::donnager::gravity::kepler as kepler;
// State
// Objects- pos, vel, mass, grav field model

pub struct GravConstants {
    pub geopotential: Geopotential,
    pub raan_dot: f64,
    pub arg_perigee_dot: f64,
    pub mean_anomaly_dot: f64,
    pub c1: f64,
    pub c4: f64,
    pub k0: f64,
    pub k1: f64,
    pub orbit_0: kepler::Orbit,
}


pub struct Geopotential {
    pub eq_radius: f64, // Equatorial radius
    pub ke: f64, // square root of earth's grav param in earth radii^3 min^-2
    pub j2: f64, // un-normalised second zonal harmonic
    pub j3: f64, // un-normalised third zonal harmonic
    pub j4: f64, // un-normalised fourth zonal harmonic
}

pub const WGS84: Geopotential = Geopotential {
    eq_radius: 6378.137,
    ke: 0.07436685316871385,
    j2: 0.00108262998905,
    j3: -0.00000253215306,
    j4: -0.00000161098761,
};


impl GravConstants {

    /// Propogate usgin sgp4 algorithm
    pub fn propogate(
        &self, 
        eval_datetimes: Vec<DateTime<Utc>>
    ) ->  Vec<Vector3<f64>> {
        let mut state: Vec<Vector3<f64>>;
        let mut time: Vec<DateTime<Utc>>;
        let mut dt: f64;
        let mut dt_max: f64;
        let mut dt_min: f64;
        let mut dt_avg: f64;
        let mut dt_sum: i64;
        let mut dt_sum_sq: f64;
        let mut dt_sum_cub: f64;
        
        dt_sum = eval_datetimes
            .iter()
            .map(|dt| dt.timestamp())
            .sum::<i64>();
        let raan = self.raan_dot * dt_sum as f64;
        let argument_of_perigee = self.arg_perigee_dot * dt_sum as f64;
        let mean_anomaly = self.mean_anomaly_dot * dt_sum as f64;
        
        let c1 = self.c1;
        let c4 = self.c4;
        let k0 = self.k0;
        let k1 = self.k1;
        let orbit = self.orbit_0.clone();
        let orbit = orbit.propogate(raan, argument_of_perigee,
            mean_anomaly, self.c1, self.c4, self.k0, self.k1);
        
        return state
    }

}