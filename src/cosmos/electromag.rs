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
/// luminosity : `f64`
///     Luminosity in W/m^2
/// 
/// radial_distance : `f64`
///     Distance from source
/// 
/// Ouputs
/// ------
/// apparent_lum : `f64`
///     Luminosity at distance
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
/// frequency : `f64`
///     Frequency of light
/// 
/// absolute_temp : `f64`
///     Temperature of blackbody in Kelvin
/// 
/// Outputs
/// -------
/// radiance : `f64`
///     Power radiated at frequency
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
/// em_flux_vec : `Vector3<f64>`
///     Electromagnetic flux vector
/// 
/// panel_area_vec : `Vector3<f64>`
///     Solar panel area normal vector
/// 
/// efficiency : `f64`
///     Solar power efficiency
/// 
/// Outputs
/// -------
/// max_power : `f64`
///     Maximum power generated
pub fn calc_max_solar_power_gen(
    em_flux_vec: Vector3<f64>,
    panel_area_vec: Vector3<f64>,
    efficiency: f64,
) -> f64 {
    let inner_prod: f64 = panel_area_vec.dot(&em_flux_vec);
    let cos_incidence: f64 = inner_prod / (panel_area_vec.norm() * em_flux_vec.norm());
    let max_power: f64 = efficiency * panel_area_vec.norm() * em_flux_vec.norm() * cos_incidence;
    return max_power
}

/// Calculate force from electromagnetic radiation
/// 
/// Inputs
/// ------
/// em_flux_vec : `f64`
///     Electromagnetic flux vector
/// 
/// obj_reflectivity
///     0.0 -> 1 -> 2.0
///     0.0 : transparent
///     1.0 : blackbody
///     2.0 : mirror
pub fn calc_em_radiative_force(
    em_flux: f64,
    src_obj_vec: Vector3<f64>,
    obj_reflecitivity: f64,
    normal_area: f64
) -> Vector3<f64> {
    let pressure = em_flux / cst::SPEED_OF_LIGHT;
    let mag = pressure * obj_reflecitivity * normal_area;
    let force = mag * src_obj_vec/ src_obj_vec.norm();
    return force
}


#[test]
fn test_solar_power(){
    let frequency = 536.34;
    let absolute_temp = 234.34;
    let radiance = calc_spectral_radiance(frequency, absolute_temp);
}
