/*
Physics calculation
*/
use crate::constants as cst;

// Calculate schwarzchild radius of a given mass
pub fn calc_schwarzchild_radius(
    mass: f64
) -> f64 {
    let radius: f64 = 2.0 * mass * cst::GRAV_CONST / cst::SPEED_OF_LIGHT.powi(2);
    return radius
}

// Calculate time dilation of relative velocity
pub fn calc_time_dilation(
    t_1: f64,
    rel_vel: f64
) -> f64 {
    let lorentz: f64 = (1.0 - rel_vel.powi(2)/cst::SPEED_OF_LIGHT.powi(2)).sqrt();
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
