/*
Physics calculation
*/
use crate::constants;

// Calculate radial ditance using two way speed of light
pub fn calc_radial_distance(
    time_delay: f64
) -> f64 {
    let r: f64 = constants::SPEED_OF_LIGHT * time_delay / 2.0;
    return r
}

// Calculate radial velocity using Doppler effect
pub fn calc_radial_vel(
    tx_wavelength: f64,
    rx_wavelength: f64
) -> f64 {
    let wavelength_shift: f64 = rx_wavelength - tx_wavelength;
    let v_r: f64 = wavelength_shift * constants::SPEED_OF_LIGHT / tx_wavelength;
    return v_r
}

// Calculate schwarzchild radius of a given mass
pub fn calc_schwarzchild_radius(
    mass: f64
) -> f64 {
    let radius: f64 = 2.0 * mass * constants::GRAV_CONST / constants::SPEED_OF_LIGHT.powi(2);
    return radius
}

// Calculate time dilation of relative velocity
pub fn calc_time_dilation(
    t_1: f64,
    rel_vel: f64
) -> f64 {
    let lorentz: f64 = (1.0 - rel_vel.powi(2)/constants::SPEED_OF_LIGHT.powi(2)).sqrt();
    let t_2: f64 = t_1 / (lorentz);
    return t_2
}

// Calculate apparant angular size of object in fov
pub fn calc_angular_size(
    object_radius: f64,
    radial_distance: f64
) -> f64 {
    let arc_rads: f64 = (object_radius / radial_distance).atan(); 
    return arc_rads
}
