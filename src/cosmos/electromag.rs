/*
Effects and interactions of the electromagnetic field
*/

use std::f64::consts::PI;
use nalgebra::Vector3;

use crate::constants as cst;

/// Calculate apparent brightness of object
/// 
/// Inputs
/// ------
/// 
pub fn calc_apparent_brightness(
    luminosity: f64,
    radial_distance: f64
) -> f64 {
    let apparent_lum: f64 = luminosity / (4.0 * PI * radial_distance.powi(2));
    return apparent_lum
}

/// Calculate spectral radiance
/// 
/// Inputs
/// ------
pub fn calc_spectral_radiance(
    frequency: f64,
    absolute_temp: f64
) -> f64 {
    let numerator: f64 = 2.0 * cst::PLANCKS_CONST * frequency.powi(3);
    let exponent: f64 = 
        cst::PLANCKS_CONST * frequency / (cst::BOLTZMAN_CONST * absolute_temp);
    let denom: f64 = cst::SPEED_OF_LIGHT.powi(2)*(exponent.exp() - 1.0);
    let radiance: f64 = numerator / denom;
    return radiance
}

/// Calculate solar power generated
/// 
/// Inputs
/// ------
pub fn calc_max_solar_power_gen(
    flux_vec: Vector3<f64>,
    area_vec: Vector3<f64>,
    efficiency: f64,
) -> f64 {
    let cos_incidence: f64 = area_vec.transpose() * flux_vec / (area_vec.norm() * flux_vec.norm());
    let max_power: f64 = efficiency * area_vec.norm() * flux_vec.norm() * cos_incidence;
    return max_power
}

pub fn calc_normal_area(
    incidence_angle: f64,

)



#[test]
fn test_solar_power(){

}
