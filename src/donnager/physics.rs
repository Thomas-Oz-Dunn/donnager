/*
Physics calculation
*/
use std::f64::consts::PI;
#[path="./constants.rs"] mod constants;

pub fn calc_radial_distance(
    time_delay: f64
) -> f64 {
    // Two way speed of light
    let r: f64 = constants::SPEED_OF_LIGHT * time_delay / 2;
    return r
}

pub fn calc_radial_vel(
    tx_wavelength: f64,
    rx_wavelength: f64
) -> f64 {
    // Doppler effect
    let wavelength_shift: f64 = rx_wavelength - tx_wavelength;
    let v_r: f64 = wavelength_shift * constants::SPEED_OF_LIGHT / tx_wavelength;
    return v_r
}

pub fn calc_schwarzchild_radius(
    mass: f64
) -> f64 {
    let radius: f64 = 2.0 * mass * constants::GRAV_CONST / constants::SPEED_OF_LIGHT.powi(2);
    return radius
}

pub fn calc_time_dilation(
    1_d_time: f64,
    rel_vel: f64
) -> f64 {
    let lorentz: f64 = (1.0 - rel_vel.powi(2)/constants::SPEED_OF_LIGHT.powi(2)).sqrt();
    let 2_d_time: f64 = 1_d_time / (lorentz);
    return 2_d_time
}

pub fn calc_angular_size(
    object_radius: f64,
    radial_distance: f64
) -> f64 {
    let arc_rads: f64 = (object_radius / radial_distance).atan(); 
    return arc_rads
}
