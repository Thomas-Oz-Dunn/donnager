/*
Gravitational Bodies
*/

use chrono::*;
use nalgebra as na;
use std::f64::consts::PI;
use na::{Vector3, Matrix3};

pub const TOLERANCE: f64 = 1e-8;

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
        let radius: f64 = self.calc_prime_vertical(lla[0]);
        let x: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].cos();
        let y: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].sin();
        let z: f64 = ((1.0 - self.eccentricity.powi(2)) * radius + lla[2]) * lla[0].sin();
        let xyz: Vector3<f64> = Vector3::new(x, y, z); 
        return xyz
    }

    // Calculate prime vertical radius to surface at latitude
    pub fn calc_prime_vertical(&self, lat_deg: f64) -> f64 {
        let lat_radians: f64 = PI * lat_deg / 180.0;
        let radius: f64 = 
            self.eq_radius / (1.0 - (self.eccentricity * lat_radians.sin()).powi(2)).sqrt();
        return radius
    }

    // Rectangular coordinates to geodetic
    // E.g. ECEF to LLH
    pub fn xyz_to_geodetic(&self, xyz: Vector3<f64>) -> Vector3<f64> {
        // Zhu's method
        let a: f64 = self.eq_radius;
        let ecc_2: f64 = self.eccentricity.powi(2);
        
        let b: f64 = (a.powi(2)*(1.0 - ecc_2)).sqrt();
        let ecc_2_prime: f64 = a.powi(2) / b.powi(2) - 1.0;
        let p: f64 = (xyz[0].powi(2) + xyz[1].powi(2)).sqrt();
        let g: f64 = p.powi(2) + (1.0 - ecc_2) * xyz[2].powi(2) - 
            ecc_2 * (a.powi(2) - b.powi(2));
        let c: f64 = ecc_2.powi(2) * 54.0 * b.powi(2) * xyz[2].powi(2) * p.powi(2) / (g.powi(3));
        let s: f64 = (1.0 + c + (c.powi(2) + 2.0 * c).sqrt()).powf(1.0 / 3.0);
        let capP: f64 = 54.0 * b.powi(2)*xyz[2].powi(2) / 
            (3.0 * (s + 1.0 + 1.0 / s).powi(2) * g.powi(2));
        let q: f64 = (1.0 + 2.0 * ecc_2.powi(2) * capP).sqrt();
        let r_0: f64 = -capP * ecc_2 * p /(1.0+q) + 
            ((a.powi(2)/2.0)*(1.0 + 1.0 / q) - 
            capP * (1.0 - ecc_2) * xyz[2].powi(2) / (q * (1.0 + q)) - 
            capP*p.powi(2)/2.0).sqrt();
        let u: f64 = ((p - ecc_2*r_0).powi(2) + xyz[2].powi(2)).sqrt();
        let v: f64 = ((p - ecc_2*r_0).powi(2) + (1.0 - ecc_2)*xyz[2].powi(2)).sqrt();
        let z_0: f64 = b.powi(2) * xyz[2] / (a * v);

        let alt: f64 = u * (1.0 - b.powi(2) / (a * v));
        let lat: f64 = ((xyz[2] + ecc_2_prime*z_0)/p).atan();
        let lon: f64 = (xyz[1] / xyz[0]).atan();
        let lla: Vector3<f64> = Vector3::new(lat, lon, alt);
        return lla
    }

}


impl SurfacePoint{

    pub fn find_time_zone(&self) -> TimeZone<Offset = > {
        let lla: Vector3<f64> = self.pos_lla;

        // Look up table?? Ideal timezones at long

        // Real timezones ughhhh
    }

    // Map between fixed frame observation to enu
    pub fn ecef_to_enu(&self, ecef: Vector3<f64>) -> Vector3<f64> {
        let pos_lla: Vector3<f64> = self.pos_lla;
        let pos_ecef: Vector3<f64> = self.Body.geodetic_to_xyz(pos_lla);
        let vec_ecef: Vector3<f64> = ecef - pos_ecef;
        let ecef_enu: Matrix3<f64> = Matrix3::new(
            -pos_lla[1].sin(), pos_lla[1].cos(), 0,
            -pos_lla[1].cos()*pos_lla[0].sin(), -pos_lla[1].sin()*pos_lla[0].sin(), pos_lla[0].cos(),
            pos_lla[1].cos()*pos_lla[0].cos(), pos_lla[1].sin()*pos_lla[0].cos(), pos_lla[0].sin());
        let enu: Vector3<f64> = ecef_enu * vec_ecef;
        return enu
    }

    // Map between enu and fixed frame
    pub fn enu_to_ecef(&self, enu: Vector3<f64>) -> Vector3<f64> {
        let enu_ecef: Matrix3<f64> = Matrix3::new(
            -pos_lla[1].sin(), -pos_lla[1].cos()*pos_lla[0].sin(), pos_lla[1].cos()*pos_lla[0].cos(),
            pos_lla[1].cos(), -pos_lla[1].sin()*pos_lla[0].sin(), pos_lla[1].sin()*pos_lla[0].cos(),
            0, pos_lla[0].cos(), pos_lla[0].sin()
        );
        let vec_ecef: Vector3<f64> = enu_ecef * enu;
        let pos_ecef: Vector3<f64> = self.Body.geodetic_to_xyz(self.pos_lla);
        let ecef: Vector<f64> = vec_ecef - pos_ecef;
        return ecef
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