

use nalgebra as na;
use na::{Vector3, Matrix3};

use crate::gravity::kepler;
use crate::constants as cst;

#[derive(Clone, Debug, PartialEq)]
pub struct SurfacePoint{
    pub name: String,
    pub body: kepler::Body,
    pub pos_lla: Vector3<f64>,
}


#[derive(Clone, Debug, PartialEq)]
pub enum Frame{
    Topocentric,            // ENU
    GeocentricFixed,        // ECEF
    GeocentricInertial,     // ECI
    Geodetic,               // LLA (Earth)
    AreocentricFixed,       // MCMF
    AreocentricInertial,    // MCI
    Areodetic,              // LLA (Mars)
    Heliocentric            // Sun
}

/// SurfacePoint methods
impl SurfacePoint{
    /// Calculate tangential velocity magnitude at surface point
    /// 
    /// Outputs
    /// -------
    /// tan_vel : `f64`
    ///     Tangential velocity magnitude in km / s
    pub fn calc_surface_vel(&self) -> f64 {
        let radius_km: f64 = self.body.eq_radius / cst::KILO;
        let equatorial_vel: f64 = self.body.rotation_rate * radius_km;
        let tan_vel: f64 = (self.pos_lla[0].cos() * equatorial_vel).abs();
        return tan_vel
    }
    
    /// Calculate radius from center of body at surface point
    /// 
    /// Outputs
    /// -------
    /// radius : `Vector3<f64>`
    ///     Radial vector in meters
    pub fn calc_surface_radius(&self) -> Vector3<f64> {
        let prime_vertical: f64 = self.body.calc_prime_vertical(self.pos_lla[0]);
        let dir: Vector3<f64> = self.body.geodetic_to_xyz(self.pos_lla);
        let radius: Vector3<f64> = (prime_vertical + self.pos_lla[2]) * dir;
        return radius
    }

    /// Calculate delta v required to reach an altitude from surface
    /// 
    /// Inputs
    /// ------
    /// altitude : `f64`
    ///     Altitude off of surface
    pub fn calc_delta_v(
        &self,
        altitude: f64,
    ) -> f64 {
        let radius: f64 = altitude + self.calc_surface_radius().norm();
        let surface_vel: f64 = self.calc_surface_vel();
        let delta_v_req: f64 = self.body.calc_orbital_velocity_mag(radius); 
        let delta_v: f64 = delta_v_req - surface_vel;
        return delta_v
    }


    /// Map between fixed frame observation to enu
    pub fn ecef_to_enu(&self, ecef_2: Vector3<f64>) -> Vector3<f64> {
        let pos_lla: Vector3<f64> = self.pos_lla;
        let pos_ecef: Vector3<f64> = self.body.geodetic_to_xyz(pos_lla);
        let vec_ecef: Vector3<f64> = ecef_2 - pos_ecef;
        let ecef_enu: Matrix3<f64> = Matrix3::new(
            -pos_lla[1].sin(), pos_lla[1].cos(), 0.0,
            -pos_lla[1].cos()*pos_lla[0].sin(), -pos_lla[1].sin()*pos_lla[0].sin(), pos_lla[0].cos(),
            pos_lla[1].cos()*pos_lla[0].cos(), pos_lla[1].sin()*pos_lla[0].cos(), pos_lla[0].sin());
        let enu: Vector3<f64> = ecef_enu * vec_ecef;
        return enu
    }

    /// Map between enu and fixed frame
    /// 
    /// 
    pub fn enu_to_ecef(&self, enu: Vector3<f64>) -> Vector3<f64> {
        let pos_lla: Vector3<f64> = self.pos_lla;
        let enu_ecef: Matrix3<f64> = Matrix3::new(
            -pos_lla[1].sin(), -pos_lla[1].cos()*pos_lla[0].sin(), pos_lla[1].cos()*pos_lla[0].cos(),
            pos_lla[1].cos(), -pos_lla[1].sin()*pos_lla[0].sin(), pos_lla[1].sin()*pos_lla[0].cos(),
            0.0, pos_lla[0].cos(), pos_lla[0].sin()
        );
        let vec_ecef: Vector3<f64> = enu_ecef * enu;
        let pos_ecef: Vector3<f64> = self.body.geodetic_to_xyz(self.pos_lla);
        let ecef: Vector3<f64> = vec_ecef - pos_ecef;
        return ecef
    }

}