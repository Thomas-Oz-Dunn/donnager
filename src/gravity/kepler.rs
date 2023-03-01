/*
Gravitational Bodies
*/

use nalgebra as na;
use chrono::DateTime as DateTime;
use std::f64::consts::PI;
use na::Vector3;

use crate::constants as cst;
use crate::gravity::barneshut as bh;

/// Gravitational Body
#[derive(Clone, Debug, PartialEq)]
pub struct Body{
    pub name: String,
    pub grav_param: f64,
    pub eq_radius: f64,
    pub rotation_rate: f64,
    pub eccentricity: f64
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


/// Orbit structure
/// 
/// Variables
/// ---------
/// name : `String`
#[derive(Clone, Debug, PartialEq)]
pub struct Orbit{
    pub name: String,
    pub grav_param: f64,
    pub semi_major_axis: f64, 
    pub eccentricity: f64,
    pub inclination: f64,
    pub raan: f64,
    pub argument_of_perigee: f64,
    pub mean_anomaly: f64,
    pub mean_motion: f64
}

impl Orbit {

    /// Populate Orbit from Keplerian parameters
    /// 
    /// Inputs
    /// ------
    /// name : `String`
    pub fn from_keplerian(
        name: String,
        grav_param: f64,
        semi_major_axis: f64, 
        eccentricity: f64,
        inclination: f64,
        raan: f64,
        argument_of_perigee: f64,
        mean_anomaly: f64,
        mean_motion: f64
    ) -> Self {
        Orbit {
            name,
            grav_param,
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            argument_of_perigee,
            mean_anomaly,
            mean_motion
        }
    }

    /// Populate Orbit from standard Two Line Element
    /// 
    /// Inputs
    /// ------
    /// tle_str : `String` 
    ///     NORAD Two Line Element Identification String
    pub fn from_tle(
        tle_str: String
    ) -> Self {
        // Two Line element usage assumes Earth Centered
        let grav_param: f64 = cst::EARTH_GRAV_PARAM;
    
        let lines: Vec<&str> = tle_str.lines().collect();
        let name: &str = lines[0];
        // let line1: Vec<&str> = lines[1].to_string().split_whitespace().collect();
        let binding: String = lines[2].to_string();
        let line2: Vec<&str> = binding.split_whitespace().collect();
        
        // let element_num: &str = line1[line1.len()];
        // let epoch: &str = line1[3];
        // let mean_motion_prime: &str = line1[4];
        // let mean_motion_2: &str = line1[5];
        let inc_str = line2[2].to_string();
        let inc: f64 = inc_str.parse::<f64>().unwrap();
        let raan: f64 = line2[3].to_string().parse::<f64>().unwrap();
        let ecc: f64 = line2[4].to_string().parse::<f64>().unwrap() * 10e-7;
        let arg_perigee: f64 = line2[5].to_string().parse::<f64>().unwrap();
        let mean_anomaly: f64 = line2[6].to_string().parse::<f64>().unwrap();
    
        let end_str: &str = line2[line2.len()-1];
        let mean_motion: f64 = end_str[..11].to_string().parse::<f64>().unwrap();
        // let rev_num: f64 = end_str[12..].to_string().parse::<f64>().unwrap();
        let semi_major_axis: f64 = (mean_motion.powi(2) / (grav_param)).powf(1.0/3.0);
    
        Orbit {
            name: name.to_string(),
            grav_param,
            semi_major_axis,
            raan,
            eccentricity: ecc,
            inclination: inc,
            argument_of_perigee: arg_perigee,
            mean_anomaly,
            mean_motion,
        }
    
    }

    /// Populate Orbit from position and velocity vectors
    /// 
    /// Inputs
    /// ------
    /// name : `String`
    ///     Name of object
    /// 
    /// grav_param : `f64`
    ///     Central body gravitational parameter
    /// 
    /// pos : `Vector3<f64>`
    ///     Position vector of object
    /// 
    /// vel : `Vector3<f64>`
    ///     Velocity vector of object
    pub fn from_pos_vel(
        name: String,
        grav_param: f64,
        pos: Vector3<f64>,
        vel: Vector3<f64>
    ) -> Self {
        let spec_ang_moment: Vector3<f64> = pos.cross(&vel);
        let spec_lin_moment: f64 = pos.dot(&vel);

        let ecc_vec: Vector3<f64> = 
            ((vel.norm().powi(2) - grav_param / pos.norm())*pos 
            - (spec_lin_moment*vel)) / grav_param;
        let node_vec: Vector3<f64> = Vector3::z_axis().cross(&spec_ang_moment);

        let semi_major_axis: f64 = 
            spec_ang_moment.norm().powi(2) * 
            (1.0 - ecc_vec.norm_squared()) / grav_param;

        Orbit {
            name,
            grav_param,
            semi_major_axis,
            eccentricity: ecc_vec.norm(),
            raan: (node_vec[0] / node_vec.norm()).acos(),
            inclination: (spec_ang_moment[2] / spec_ang_moment.norm()).acos(),
            argument_of_perigee: node_vec.angle(&ecc_vec),
            mean_anomaly: ecc_vec.angle(&pos),
            mean_motion: 1.0 / (2.0 * PI * (semi_major_axis.powi(3)/grav_param).sqrt())
        }

    }

    /// Calculate object position
    /// 
    /// Inputs
    /// ------
    /// eval_datetime
    /// 
    /// method
    pub fn calc_position(
        self, 
        eval_datetime: DateTime,
        method: String
    ) -> Vec<f64> {

        let valid_methods = [
            "runge_kutta", "kepler", "barnes_hut"];
        assert!(valid_methods.contains(&method));
        let mut pos;

        if method == "kepler"{

            // M = E - e sin E


            pos = vec![0.,0.,0.];
        } else if method == "barnes_hut"{

            let motion_0 = ;
            let satellite: Particle = Particle { mass: (), motion: motion_0};


            let earth: Body = Body {
                name: "Earth".to_string(),
                grav_param: cst::EARTH_GRAV_PARAM,
                eq_radius: cst::EARTH_RADIUS_EQUATOR,
                rotation_rate: cst::EARTH_ROT_RATE,
                eccentricity: cst::EARTH_ECC
            };
 
            let earth_motion: Vec<Vector3<f64>> = vec![Vector3::zeros(); 3];
            let earth_particle: Particle = earth.to_particle(earth_motion);
        
            let mut particles: Box<[Particle]> = vec![earth_particle, satellite].into_boxed_slice();
    
            let step_size: f64 = 0.1;
            let theta: f64 = 1.0;
            let n_steps: usize = 10;
            let is_debug: bool = false;
        
            particles = bh::barnes_hut_gravity(
                particles, step_size, n_steps, theta, is_debug);
            pos = particles[1].motion[0];
        
        }



        return pos
    }

    /// Next overhead pass
    /// 
    /// Inputs
    /// ------
    /// pos_lla
    /// 
    /// 
    pub fn calc_next_overpass(
        self, 
        pos_lla: Vector3<f64>, 
        overhead_cone_angle: f64
    ) -> Datetime {

        // use period and inclination to get lattitude sine
        // use period and earth rotation to get long sine
    }

}


/// Calculate total delta v for hohmann transfer
/// 
/// Inputs
/// ------
/// radius_1 : `f64`
///     Radius of inner orbit
/// 
/// radius_2 : `f64`
///     Radius of outer orbit
/// 
/// vel_0 : `f64`
///     Initial Velocity of vehicle
/// 
/// Outputs
/// -------
/// delta_v_total : `f64`
///     Total delta v for maneuver
pub fn calc_hohmann_transfer(
    radius_1: f64,
    radius_2: f64,
    vel_0: f64
) -> f64 {
    let delta_v_1: f64 = vel_0 * (((2.0 * radius_2)/(radius_1 + radius_2)).sqrt()- 1.0);
    let delta_v_2_num: f64 = (radius_1 / radius_2).sqrt() + (2.0 * radius_1);
    let delta_v_2_den: f64 = (radius_2 *(1.0 + radius_2 / radius_1)).sqrt();
    let delta_v_2: f64 = vel_0 * delta_v_2_num / delta_v_2_den; 
    let delta_v_total: f64 = delta_v_1 + delta_v_2;
    return delta_v_total
}


/// Calculate sphere of influence for a body
/// 
/// Inputs
/// ------
/// mass_0 : `f64`
///     Primary mass of 2 body system
/// 
/// mass_1 : `f64`
///     Secondary mass of 2 body system
/// 
/// semi_major_axis : `f64`
///     Orbital ellipse semi major axis
/// 
/// eccentricity : `f64`
///     Orbital ellipse eccentricity
/// 
/// Outputs
/// -------
/// radius : `f64`
///     Radius of sphere of influence
pub fn calc_hill_sphere(
    mass_0: f64,
    mass_1: f64,
    semi_major_axis: f64,
    eccentricity: f64
) -> f64 {
    let radius: f64 = 
        semi_major_axis * (1.0 - eccentricity) * (mass_0 / (3.0 * mass_1)).powf(1.0/3.0);
    return radius
}


#[cfg(test)]
mod orbit_tests {
    use super::*;

    #[test]
    fn test_hohmann_transfer(){
        let radius_1: f64 = 5000.;
        let radius_2: f64 = 5500.;
        let vel_0: f64 = 10.;
        let deltav: f64 = calc_hohmann_transfer(radius_1, radius_2, vel_0);

        assert_eq!(deltav, 930.8082549013038);

    }

    #[test]
    fn test_hill_sphere(){
        let earth_mass: f64 = crate::constants::EARTH_MASS;
        let sun_mass: f64 = crate::constants::SUN_MASS;
        let earth_orbit_semi_major: f64 = crate::constants::EARTH_ORBIT_SEMI_MAJOR;
        let earth_orbit_ecc: f64 = crate::constants::EARTH_ORBIT_ECC;

        let sphere_rad: f64 = calc_hill_sphere(
            sun_mass, 
            earth_mass,
            earth_orbit_semi_major,
            earth_orbit_ecc
        );

        assert_eq!(sphere_rad, 7069993624.788241)
    }

    #[test]
    fn test_orbit(){
        let tle_str = "ISS (ZARYA) 
        1 25544U 98067A   23035.69666365  .00008902  00000+0  16600-3 0  9994
        2 25544  51.6420 264.7747 0008620 314.4274 150.8239 15.49588766381243";
        let orb = Orbit::from_tle(tle_str.to_string());

        assert_eq!(orb.argument_of_perigee, 314.4274);
    }


}
