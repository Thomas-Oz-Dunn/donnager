

use std::f64::consts::PI;
use crate::constants as cst;

pub fn calc_apparent_brightness(
    luminosity: f64,
    radial_distance: f64
) -> f64 {
    let apparent_lum: f64 = luminosity / (4.0 * PI * radial_distance.powi(2));
    return apparent_lum
}

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
