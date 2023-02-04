/*
Astrodynamics
*/

use nalgebra as na;
use std::f64::consts::PI;
use na::Vector3;

use crate::constants;

/// Orbit structure
/// 
/// Variables
/// ---------
/// 
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
        let grav_param: f64 = constants::EARTH_GRAV_PARAM;
    
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
            semi_major_axis: semi_major_axis,
            raan: raan,
            eccentricity: ecc,
            inclination: inc,
            argument_of_perigee: arg_perigee,
            mean_anomaly: mean_anomaly,
            mean_motion: mean_motion,
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
            name: name,
            grav_param: grav_param,
            semi_major_axis: semi_major_axis,
            eccentricity: ecc_vec.norm(),
            raan: (node_vec[0] / node_vec.norm()).acos(),
            inclination: (spec_ang_moment[2] / spec_ang_moment.norm()).acos(),
            argument_of_perigee: node_vec.angle(&ecc_vec),
            mean_anomaly: ecc_vec.angle(&pos),
            mean_motion: 1.0 / (2.0 * PI * (semi_major_axis.powi(3)/grav_param).sqrt())
        }

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
        let radius_1 = 5000.;
        let radius_2 = 5500.;
        let vel_0 = 10.;
        let deltav = calc_hohmann_transfer(radius_1, radius_2, vel_0);

        assert_eq!(deltav, 930.8082549013038);

    }

    #[test]
    fn test_hill_sphere(){
        let earth_mass = crate::constants::EARTH_MASS;
        let sun_mass = crate::constants::SUN_MASS;
        let earth_orbit_semi_major = crate::constants::EARTH_ORBIT_SEMI_MAJOR;
        let earth_orbit_ecc = crate::constants::EARTH_ORBIT_ECC;

        let sphere_rad = calc_hill_sphere(
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
