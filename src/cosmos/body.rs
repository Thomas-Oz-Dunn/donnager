/*
Gravitational Bodies
*/

use chrono::*;
use nalgebra as na;
use std::f64::consts::PI;
use na::{Vector3, Matrix3};


#[derive(Clone, Debug, PartialEq)]
pub struct Body{
    pub name: String,
    pub grav_param: f64,
    pub eq_radius: f64,
    pub rotation_rate: f64,
}

pub struct SurfacePoint{
    pub Body: Body,
    pub pos_lla: Vector3<f64>,
}

impl Body {

    // Calculate gravitational acceleration at radial distance
    pub fn calc_grav_acc(
        &self, 
        radius: f64
    ) -> f64 {
        let grav_acc: f64 = self.grav_param / radius.powi(2);
        return grav_acc
    }
    

    // Calculate required orbital velocity at radial distance
    pub fn calc_orbital_velocity(
        &self,
        radius: f64
    ) -> f64 {
        // TODO-TD: Vectorize
        let vel: f64 = (2.0 * self.grav_param / radius).sqrt();
        return vel
    }

    // Calculate period of orbit
    pub fn calc_period(
        &self,
        semi_major_axis: f64
    ) -> f64 {
        let time: f64 = 2.0 * PI * (semi_major_axis.powi(3)/self.grav_param).sqrt();
        return time
    }
    
    // Calculate radius for stationary orbit above body surface
    pub fn calc_stationary_orbit(&self) -> f64 {
        let period: f64 = 2.0 * PI / self.rotation_rate;
        let a: f64 = self.grav_param * period.powi(2); // 
        let r_mag: f64 = (a / (4.0 * PI.powi(2))).powf(1.0 / 3.0); // 
        return r_mag
    }

    // Calculate tangential velocity on surface of body
    pub fn calc_surface_vel(
        &self,
        pos_llh: Vector3<f64>
    ) -> f64 {
        let equatorial_vel: f64 = self.rotation_rate * self.eq_radius;
        let tan_vel: f64 = (pos_llh[0].cos() * equatorial_vel).abs();
        return tan_vel
    }

    // Calculate transformation matrix from fixed to inertial frame
    // E.g. ECEF to ECI
    pub fn calc_fixed_to_inertial_matrix(
        &self,
        datetime_utc: NaiveDateTime
    ) -> Matrix3<f64> {
        // Datetime to julian time

        // Sidereal time

        // Turn angle
        
        // Assemble matrix from angle
    }

    pub fn calc_lla_to_fixed_matrix(&self) -> Matrix3<f64> {
        // Ellipsoid, how detailed?
    }
}


impl SurfacePoint{

    pub fn find_time_zone(&self) -> TimeZone<Offset = > {
        let lla: Vector3<f64> = self.pos_lla;

        // Look up table?? Ideal timezones at long

        // Real timezones ughhhh
    }

    // Mapping between fixed frame and enu
    pub fn enu_matrix(&self, datetime_utc: NaiveDateTime) {
        let lla: Vector3<f64> = self.pos_lla;

        let fixed: Vector3<f64> = ;

        // Up == radial direction

        // North == 90 degree from lattitude

        // East == 90 degrees from longitude
    }

    pub fn next_overhead_pass(
        &self, 
        orbit: Orbit, 
        datetime_utc: NaiveDateTime
    ) -> NaiveDateTime {
        let lla: Vector3<f64> = self.pos_lla;

        // Current datetime to julian

        // Pos vel at current julian

        // ENU at current julian

        // Proportion step to angle between Up vector and pos vector


    }
}