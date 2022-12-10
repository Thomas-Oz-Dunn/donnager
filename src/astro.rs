/*
Astrodynamics
*/

use nalgebra as na;
use std::f64::consts::PI;
use na::Vector3;

#[path="./constants.rs"] mod constants;

pub struct Orbit{
    name: String,
    semi_major_axis: f64, 
    eccentricity: f64,
    inclination: f64,
    argument_of_perigee: f64,
    mean_anomaly: f64,
    mean_motion: f64
}

pub fn read_tle(
    tle_str: String
) -> Orbit {
    // parse the string
    // Example
    // ISS (ZARYA)
    // 1 25544U 98067A   20045.18587073  .00000950  00000-0  25302-4 0  9990
    // 2 25544  51.6443 242.0161 0004885 264.6060 207.3845 15.49165514212791
    // name: String::from("ISS (ZARYA)"),
    // satellite_number: 25544,
    // classification: 'U',
    // international_designator: String::from("98067A"),
    // epoch: String::from("20045.18587073"),
    // first_derivative_mean_motion: 0.00000950,
    // second_derivative_mean_motion: 0.0,
    // drag_term: 0.25302e-4,
    // ephemeris_type: 0,
    // element_number: 999,
    // inclination: 51.6443,
    // right_ascension: 242.0161,
    // eccentricity: 0.0004885,
    // argument_of_perigee: 264.6060,
    // mean_anomaly: 207.3845,
    // mean_motion: 15.49165514,
    // revolution_number: 21279,

    let tle_orbit = Orbit {
        name: name,
        semi_major_axis: a,
        eccentricity: e,
        inclination: i,
        argument_of_perigee: o,
        mean_anomaly: u,
        mean_motion: v,
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
    let r: f64 = (grav_param.sqrt() * period / (2.0 * PI)).powf(2.0 / 3.0);
    return r
}

pub fn calc_sphere_of_influence(
    mass_0: f64,
    mass_1: f64,
    distance: f64
) -> f64 {
    let radius: f64 = distance / (1.0 + (mass_1 / mass_0).sqrt());
    return radius
}

pub fn calc_orbit_parameters(
    name: String,
    grav_param: f64,
    pos: Vector3<f64>,
    vel: Vector3<f64>
) -> Orbit {
    // Turn into singular matrix transformation?
    // Wedge product?
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
        inclination: (spec_ang_moment[2] / spec_ang_moment.norm()).acos(),
        argument_of_perigee: node_vec.angle(&ecc_vec),
        mean_anomaly: ecc_vec.angle(&pos),
        mean_motion: 1.0 / calc_period(grav_param, semi_major_axis)
    };

    return orbit
}
