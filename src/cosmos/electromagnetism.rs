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

/// Calculate wavelength of peak black body radiance (m)
/// 
/// Inputs
/// ------
/// absolute_temp : `f64`
///     Absolute blackbody temperature (K)
/// 
/// Outputs
/// -------
/// wavelength : `f64`
///     Wavelength of peak radiance
pub fn calc_peak_wavelength(
    absolute_temp: f64
) -> f64 {
    let wiens_dis_const = 2.897771955;
    let wavelength = wiens_dis_const / absolute_temp;
    return wavelength
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


/// Calculate radial ditance using two way speed of light
/// 
/// Inputs
/// ------
/// time_delay : `f64`
///     Time delay in seconds
/// 
/// Outputs
/// -------
/// c : `f64`
///     Radial magnitudes
pub fn calc_radial_distance(
    time_delay: f64
) -> f64 {
    let radius: f64 = cst::SPEED_OF_LIGHT * time_delay / 2.0;
    return radius
}

/// Calculate radial velocity using Doppler effect
/// 
/// Inputs
/// ------
/// 
pub fn calc_radial_vel(
    tx_wavelength: f64,
    rx_wavelength: f64
) -> f64 {
    let wavelength_shift: f64 = rx_wavelength - tx_wavelength;
    let v_r: f64 = wavelength_shift * cst::SPEED_OF_LIGHT / tx_wavelength;
    return v_r
}


#[cfg(test)]
mod electromag_tests{
    use super::*;

    #[test]
    fn test_radiance(){
        let absolute_temp = 234.34;
        let frequency = cst::SPEED_OF_LIGHT / calc_peak_wavelength(absolute_temp);
        let radiance = calc_spectral_radiance(frequency, absolute_temp);
        assert_eq!(radiance, 4.221293395184273e-17)
    
    }

    #[test]
    fn test_power(){
        let efficiency = 0.1;
        let solar_flux_mag = cst::SUN_MEAN_SOLAR_FLUX;
        let em_flux_vec: Vector3<f64> = Vector3::new(0., 0., 1.) * solar_flux_mag;
        let panel_area_vec: Vector3<f64> = 
            Vector3::new(0., 0., 1.);
        let max_power = calc_max_solar_power_gen(
            em_flux_vec,
            panel_area_vec,
            efficiency
        );
        assert_eq!(max_power, 136.691)
    }
}

