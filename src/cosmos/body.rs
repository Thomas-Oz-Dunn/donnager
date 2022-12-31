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
    pub eccentricity: f64
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

    // Transform from fixed to inertial frame
    // E.g. ECEF to ECI
    pub fn fixed_to_inertial(
        &self,
        ecef: Vector3<f64>,
        datetime_utc: NaiveDateTime
    ) -> Vector3<f64> {
        // Datetime to julian date
        datetime_utc.
        // Days since 2000 January 1, 12 h UTl?
        // GMST
        // Local Sidereal time = Greenwich sidereal time + lattitude 
        let lcl_sidereal_time: f64 = gmst + ecef;

        // Precession

        // Nutation 

        // Turn angle
        
    }

    // Geodetic to rectangular coordinates
    // E.g. Latitude, Longitude, Altitude to ECEF
    pub fn geodetic_to_xyz(&self, lla: Vector3<f64>) -> Vector3<f64> {
        let semi_major: f64 = self.eq_radius;
        let eccentricity: f64 = self.eccentricity;
        let radius: f64 = semi_major / (1.0 - (eccentricity * lla[0].sin()).powi(2)).sqrt();
        
        let x: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].cos();
        let y: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].sin();
        let z: f64 = ((1.0 - eccentricity.powi(2)) * radius + lla[2]) * lla[0].sin();
        let xyz: Vector3<f64> = Vector3::new(x, y, z); 
        return xyz
    }

    // Rectangular coordinates to geodetic
    // E.g. ECEF to LLH
    pub fn xyz_to_geodetic(&self, xyz: Vector3<f64>) -> Vector3<f64> {
        let semi_major: f64 = self.eq_radius;
        let ecc: f64 = self.eccentricity;

        let radius: f64 = xyz.magnitude();

        // x = a * cos(lat) * cos(lon) / sqrt(1 - e^2sin^2(lat)) + alt)
        // y = a * cos(lat) * sin(lon) / sqrt(1 - e^2sin^2(lat)) + alt)
        // z = ((1 - e^2) * R + alt) * sin(lat)
        // R = a / sqrt(1 - e^2*sin^2(lat))
        // z = ((1 - e^2) * a / (sqrt(1 - e^2*sin^2(lat))) + alt) * sin(lat)

        // x / y = cos(lon) / sin(lon) 
        // y / x = tan(lon)
        // lon = atan(y/x)


        // z / y = sin(lat) / cos(lat) * ((1 - e^2) * a / (sqrt(1 - e^2*sin^2(lat))) + alt) / (a * sin(lon) / sqrt(1 - e^2sin^2(lat)) + alt))
        // z / y = tan(lat) * ((1-e^2) * R + alt)) / 
                        //      (sin(atan(y/x)) * (radius))
        // z / y = tan(lat) * (1 - e^2 + alt / R) / (sin(atan(y/x)))
        // lat = atan(z/y * (sin(atan(y/x)) / (1 - e^2 + alt / R))

        let prime_verical:f64 = ;

        let x: f64 = (xyz[1] / xyz[0]).atan();
        let longitude: f64 = ;
        let altitude: f64 = radius - prime_vertical;
        let lla: Vector3<f64> = Vector3::new(x, longitude, altitude);

        return lla
    }

}


impl SurfacePoint{

    pub fn find_time_zone(&self) -> TimeZone<Offset = > {
        let lla: Vector3<f64> = self.pos_lla;

        // Look up table?? Ideal timezones at long

        // Real timezones ughhhh
    }

    // Mapping between fixed frame and enu
    pub fn enu_ecef(&self, datetime_utc: NaiveDateTime) {
        let pos_lla: Vector3<f64> = self.pos_lla;
        let pos_ecef: Vector3<f64> = self.Body.geodetic_to_rect(pos_lla);

        // Assumes spherical, handle z?
        let up: Vector3<f64> = pos_ecef / pos_ecef.magnitude();

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

        // Current datetime to julian date

        // Calculate orbital pos vel at current julian

        // ENU of surface point at current julian

        // Proportion step to angle between Up vector and pos vector


    }
}