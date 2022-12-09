/*
Astrodynamics Functions
*/

use std::f64::consts::PI;
use nalgebra as na;
use na::Vector3;

#[path="./constants.rs"] mod constants;

pub struct Body{
    name: String,
    mass: f64,
    radius: f64,
    period: f64
}

pub struct Orbit{
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
    let a: f64 = 4.352;
    let e: f64 = 4.352;
    let i: f64 = 4.352;
    let o: f64 = 4.352;
    let u: f64 = 4.352;
    let v: f64 = 4.352;

    let tle_orbit = Orbit {
        semi_major_axis: a,
        eccentricity: e,
        inclination: i,
        argument_of_perigee: o,
        mean_anomaly: u,
        mean_motion: v,
    };

    return tle_orbit
}


pub fn calc_grav_acc(
    mass: f64,
    radius: f64
) -> f64{
    let grav_acc: f64 = constants::GRAV_CONST * mass / radius.powi(2);
    return grav_acc
}

pub fn calc_orbital_velocity(
    mass: f64,
    radius: f64
) -> f64 {
    let vel: f64 = (2.0 * constants::GRAV_CONST * mass / radius).sqrt();
    return vel
}

pub fn calc_period(
    mass: f64,
    semi_major_axis: f64
) -> f64 {
    let grav_param: f64 = mass * constants::GRAV_CONST;
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
    mass: f64,
    period: f64
) -> f64 {
    let grav_param: f64 = mass * constants::GRAV_CONST;
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
    mass: f64,
    pos: Vector3<f64>,
    vel: Vector3<f64>
) -> Orbit {
    // Turn into singular matrix transformation?
    let grav_param: f64 = mass * constants::GRAV_CONST;

    // Wedge product?
    let spec_ang_moment: Vector3<f64> = pos.cross(&vel);
    let spec_lin_moment: f64 = pos.dot(&vel);

    let ecc_vec: Vector3<f64> = 
        ((vel.norm().powi(2) - grav_param / pos.norm())*pos 
        - (spec_lin_moment*vel)) / grav_param;
    let node_vec: Vector3<f64> = Vector3::z_axis().cross(&ang_momentum);

    let semi_major_axis: f64 = 
    spec_ang_moment.norm().powi(2) * (1.0 - ecc_vec.norm_squared()) / grav_param;

    let orbit: Orbit = Orbit {
        semi_major_axis: semi_major_axis,
        eccentricity: ecc_vec.norm(),
        inclination: (spec_ang_moment[2] / spec_ang_moment.norm()).acos(),
        argument_of_perigee: node_vec.angle(&ecc_vec),
        mean_anomaly: ecc_vec.angle(&pos),
        mean_motion: 1.0 / calc_period(mass, semi_major_axis)
    };

    return orbit
}
