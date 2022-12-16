/*
Astrodynamics
*/

use nalgebra as na;
use std::f64::consts::PI;
use na::{Vector3};

#[path="./constants.rs"] mod constants;

pub struct Orbit{
    pub name: String,
    pub semi_major_axis: f64, 
    pub eccentricity: f64,
    pub inclination: f64,
    pub raan: f64,
    pub argument_of_perigee: f64,
    pub mean_anomaly: f64,
    pub mean_motion: f64
}

pub fn read_tle(
    tle_str: String
) -> Orbit {
    // TODO- Create test
    let grav_param: f64 = constants::GRAV_CONST * constants::EARTH_MASS;

    let lines: Vec<&str> = tle_str.lines().collect();
    let name: &str = lines[0];
    let line1: Vec<&str> = lines[1].to_string().split_whitespace().collect();
    let line2: Vec<&str> = lines[2].to_string().split_whitespace().collect();

    let element_num: &str = line1[line1.len()];
    let epoch: &str = line1[3];
    let mean_motion_prime: &str = line1[4];
    let mean_motion_2: &str = line1[5];
    let drag: &str = lines[6];
    let inc: f64 = line2[2].to_string().parse::<f64>().unwrap();
    let raan: &str = line2[3].to_string().parse::<f64>().unwrap();
    let ecc: f64 = line2[4].to_string().parse::<f64>().unwrap() * 10e-7;
    let arg_perigee: &str = line2[5].to_string().parse::<f64>().unwrap();
    let mean_anomaly: &str = line2[6].to_string().parse::<f64>().unwrap();

    let end_str: &str = line2[line2.len()];
    let mean_motion: f64 = end_str[..11].to_string().parse::<f64>().unwrap();
    let rev_num: f64 = end_str[12..].to_string().parse::<f64>().unwrap();
    let semi_major_axis: f64 = (mean_motion.powi(2) / (grav_param)).powf(1.0/3.0);

    let tle_orbit: Orbit = Orbit {
        name: name.to_string(),
        semi_major_axis: semi_major_axis,
        raan: raan,
        eccentricity: ecc,
        inclination: inc,
        argument_of_perigee: arg_perigee,
        mean_anomaly: mean_anomaly,
        mean_motion: mean_motion,
    };

    return tle_orbit
}

pub fn calc_orbital_velocity(
    grav_param: f64,
    radius: f64
) -> f64 {
    let vel: f64 = (2.0 * grav_param / radius).sqrt();
    return vel
}

pub fn calc_period(
    grav_param: f64,
    semi_major_axis: f64
) -> f64 {
    let time: f64 = 2.0 * PI * (semi_major_axis.powi(3)/grav_param).sqrt();
    return time
}

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

pub fn calc_stationary_orbit(
    grav_param: f64,
    period: f64
) -> f64 {
    let r_mag: f64 = (grav_param.sqrt() * period / (2.0 * PI)).powf(2.0 / 3.0);
    return r_mag
}

pub fn calc_sphere_of_influence(
    mass_0: f64,
    mass_1: f64,
    distance: f64
) -> f64 {
    let radius: f64 = distance / (1.0 + (mass_1 / mass_0).sqrt());
    return radius
}

pub fn calc_surface_vel(
    rotation_rate: f64,
    equatorial_radius: f64,
    pos_llh: Vector3<f64>
) -> f64 {
    let equatorial_vel: f64 = rotation_rate * equatorial_radius;
    let tan_vel: f64 = pos_llh[0].cos() * equatorial_vel;
    return tan_vel
}

pub fn calc_orbit_parameters(
    name: String,
    grav_param: f64,
    pos: Vector3<f64>,
    vel: Vector3<f64>
) -> Orbit {
    let spec_ang_moment: Vector3<f64> = pos.cross(&vel);
    let spec_lin_moment: f64 = pos.dot(&vel);

    let ecc_vec: Vector3<f64> = 
        ((vel.norm().powi(2) - grav_param / pos.norm())*pos 
        - (spec_lin_moment*vel)) / grav_param;
    let node_vec: Vector3<f64> = Vector3::z_axis().cross(&spec_ang_moment);

    let semi_major_axis: f64 = 
    spec_ang_moment.norm().powi(2) * (1.0 - ecc_vec.norm_squared()) / grav_param;

    let orbit: Orbit = Orbit {
        name: name,
        semi_major_axis: semi_major_axis,
        eccentricity: ecc_vec.norm(),
        raan: (node_vec[0] / node_vec.norm()).acos(),
        inclination: (spec_ang_moment[2] / spec_ang_moment.norm()).acos(),
        argument_of_perigee: node_vec.angle(&ecc_vec),
        mean_anomaly: ecc_vec.angle(&pos),
        mean_motion: 1.0 / calc_period(grav_param, semi_major_axis)
    };

    return orbit
}