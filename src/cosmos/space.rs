

use nalgebra as na;
use na::{Vector3, Matrix3};
use std::f64::consts::PI;

use crate::constants as cst;

/// Gravitational Body
#[derive(Clone, Debug, PartialEq)]
pub struct Body{
    pub name: String,
    pub grav_param: f64,
    pub eq_radius: f64,
    pub rotation_rate: f64,
    pub eccentricity: f64
}

pub enum ReferenceFrames {
    ECI,
    ECEF,
    PFCL,
}

/// Gravitational Body
impl Body {

    /// Convert into simple particle type
    /// 
    /// Inputs
    /// ------
    /// motion: `Vec<Vector3<f64>>`
    ///     motion of body in frame
    /// 
    /// Outputs
    /// -------
    /// particle: `Particle`
    ///     `Particle` equivalent of body
    pub fn to_particle(&self, motion: Vec<Vector3<f64>>) -> Particle {
        let mass: f64 = self.grav_param / cst::GRAV_CONST;
        let particle: Particle = Particle {mass: mass, motion: motion};
        particle
    }


    /// Calculate gravitational acceleration at radial distance from Body
    /// 
    /// Inputs
    /// ------
    /// radius: `Vector3<f64>`
    ///     Radius in km from Body center
    /// 
    /// Outputs
    /// -------
    /// grav_acc: `Vector3<f64>`
    ///     Acceleration rate due to gravity
    pub fn calc_body_grav(&self, radius: Vector3<f64>) -> Vector3<f64> {
        let grav_acc: Vector3<f64>  = self.grav_param * radius / radius.norm().powi(3);
        return grav_acc
    }
    

    /// Calculate required orbital velocity at radial distance
    /// 
    /// Inputs
    /// ------
    /// radius: `f64`
    ///     Radius in km from Body center
    ///     
    /// Outputs
    /// -------
    /// val: `f64`
    ///     Required tangential velocity magnitude
    pub fn calc_orbital_velocity_mag(&self, radius: f64) -> f64 {
        let vel: f64 = (2.0 * self.grav_param / radius).sqrt();
        return vel
    }

    /// Calculate period of orbit at semi major axis
    /// 
    /// Inputs
    /// ------
    /// semi_major_axis: `f64`
    ///     Semi major axis of orbital ellipse
    /// 
    /// Outputs
    /// -------
    /// period: `f64`
    ///    Period of orbit in seconds
    pub fn calc_period(&self, semi_major_axis: f64) -> f64 {
        let period: f64 = 2.0 * PI * (semi_major_axis.powi(3)/self.grav_param).sqrt();
        return period
    }
    
    /// Calculate radius for stationary orbit above body surface
    /// 
    /// Outputs
    /// -------
    /// r_mag: `f64`
    ///     Magnitude of radius for stationary orbit
    pub fn calc_stationary_orbit(&self) -> f64 {
        let period: f64 = 2.0 * PI / self.rotation_rate;
        let a: f64 = self.grav_param * period.powi(2); 
        let r_mag: f64 = (a / (4.0 * PI.powi(2))).powf(1.0 / 3.0); 
        return r_mag
    }

    /// Geodetic to rectangular coordinates
    /// E.g. Latitude, Longitude, Altitude to ECEF
    /// 
    /// Inputs
    /// ------
    /// lla: `Vector3<f64>`
    ///     geodetic coords
    /// 
    /// Outputs
    /// -------
    /// xyz: `Vector3<f64>`
    ///     Cartesian coords
    pub fn geodetic_to_xyz(&self, lla: Vector3<f64>) -> Vector3<f64> {
        let radius: f64 = self.calc_prime_vertical(lla[0]);
        let x: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].cos();
        let y: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].sin();
        let z: f64 = ((1.0 - self.eccentricity.powi(2)) * radius + lla[2]) * lla[0].sin();
        let xyz: Vector3<f64> = Vector3::new(x, y, z); 
        return xyz
    }

    /// Calculate prime vertical radius to surface at latitude
    /// 
    /// Inputs
    /// ------
    /// lat_deg: `f64`
    ///     Lattitude in degrees
    pub fn calc_prime_vertical(&self, lat_deg: f64) -> f64 {
        let lat_radians: f64 = PI * lat_deg / 180.0;
        let radius: f64 = 
            self.eq_radius / (1.0 - (self.eccentricity * lat_radians.sin()).powi(2)).sqrt();
        return radius
    }

    /// Rectangular coordinates to geodetic
    /// E.g. ECEF to LLH
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
        let cap_p: f64 = 54.0 * b.powi(2)*xyz[2].powi(2) / 
            (3.0 * (s + 1.0 + 1.0 / s).powi(2) * g.powi(2));
        let q: f64 = (1.0 + 2.0 * ecc_2.powi(2) * cap_p).sqrt();
        let r_0: f64 = -cap_p * ecc_2 * p /(1.0+q) + 
            ((a.powi(2)/2.0)*(1.0 + 1.0 / q) - 
            cap_p * (1.0 - ecc_2) * xyz[2].powi(2) / (q * (1.0 + q)) - 
            cap_p*p.powi(2)/2.0).sqrt();
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



/// Gravitational Particle
#[derive(Clone, Debug, PartialEq)]
pub struct Particle {
    pub mass: f64,
    pub motion: Vec<Vector3<f64>>,
}

impl Particle {

    /// Convert into Body type
    /// 
    /// Inputs
    /// ------
    /// name : `String`
    ///     Name of body
    /// 
    /// eq_radius : `f64`
    ///     Equatorial radius of body
    /// 
    /// rotation_rate : `f64`
    ///     Rotation rate of body
    /// 
    /// eccentricity : `f64`
    ///     Body oblateness
    pub fn to_body(&self, name: String, eq_radius: f64, rotation_rate: f64, eccentricity: f64) -> Body {
        let grav_param: f64 = self.mass * cst::GRAV_CONST;
        let body: Body = Body {
            name,
            grav_param,
            eq_radius,
            rotation_rate,
            eccentricity
        };
        body
    }
}


#[derive(Clone, Debug, PartialEq)]
pub struct SurfacePoint{
    pub name: String,
    pub body: Body,
    pub pos_lla: Vector3<f64>,
}


#[derive(Clone, Debug, PartialEq)]
pub enum Frame{
    ENU,            // ENU
    ECEF,        // ECEF
    ECI,     // ECI
    LlaE,               // LLA (Earth)
    MCMF,       // MCMF
    MCI,    // MCI
    LlaM,              // LLA (Mars)
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