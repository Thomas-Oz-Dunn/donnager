/*
Fluid dynamics
*/

use crate::donnager_rs::constants as cst;

/// Calculate mass flow of nozzle
/// 
/// Inputs
/// ------
/// throat_area : `f64`
///     Area of nozzle throat
/// 
/// p_chamber : `f64`
///     Chamber pressure
/// 
/// c_star : `f64`
///     Characteristic velocity
/// 
/// Outputs
/// -------
/// mass_flow : `f64`
///     Mass flow rate through throat
pub fn calc_mass_flow(
    throat_area: f64,
    p_chamber: f64,
    c_star: f64
) -> f64 {
    let mass_flow: f64 = p_chamber * throat_area / c_star; 
    return mass_flow
}


pub fn calc_atmos_pressure(
    p_ref: f64,
    t_ref: f64,
    delta_h: f64,
    molar_mass: f64
) -> f64 {
    let exp: f64 = 
        cst::GRAV_CONST * molar_mass * delta_h / (cst::GAS_CONST * t_ref);
    let pressure: f64 = p_ref * exp.exp();
    return pressure
}


/// Calculate mach number
/// 
/// Inputs
/// ------
/// velocity : `f64`
///     Fluid or object velocity
/// 
/// temperature : `f64`
///     Fluid temperature
/// 
/// heat_capacity_ratio : `f64`
///     Heat capacity ratio of fluid
/// 
/// Outputs
/// -------
/// mach_number : `f64`
///     Fluid or object  mach number
pub fn calc_mach_number(
    velocity: f64,
    temperature: f64,
    heat_capacity_ratio: f64
) -> f64 {
    let speed_of_sound: f64 = 
        (heat_capacity_ratio * cst::GAS_CONST * temperature).sqrt();
    let mach_number: f64 = velocity / speed_of_sound;
    return mach_number 
}


pub fn calc_aerodynamic_force(
    density: f64,
    velocity: f64,
    coeff: f64,
    normal_area: f64
) -> f64 {
    // Can be used for drag or lift or both
    let force: f64 = density * velocity.powi(2) * coeff * normal_area / 2.0;
    return force
}


pub fn calc_characteristic_vel(
    heat_capacity_ratio: f64,
    t_chamber: f64
) -> f64 {
    let c_f: f64 = 
        (2.0/(heat_capacity_ratio + 1.0)).powf(
            -(heat_capacity_ratio + 1.0)/(2.0*(heat_capacity_ratio - 1.0))); 
    let exhaust_vel: f64 = 
        (cst::GAS_CONST * t_chamber * heat_capacity_ratio).sqrt() / 
        heat_capacity_ratio;
    let c_star: f64 = exhaust_vel * c_f;
    return c_star
}


pub fn calc_thrust_coeff(
    heat_capacity_ratio: f64, 
    expansion_ratio: f64,
    p_chamber: f64,
    p_atm: f64,
    p_exhaust: f64
) -> f64 {
    let a: f64 = 
        2.0 * heat_capacity_ratio.powi(2) / (heat_capacity_ratio - 1.0);
    let b: f64 = 
        (2.0 / (heat_capacity_ratio + 1.0)).powf(
            (heat_capacity_ratio + 1.0)/(heat_capacity_ratio - 1.0));
    let c: f64 = 
        1.0 - (p_exhaust / p_chamber).powf(
            (heat_capacity_ratio - 1.0) / heat_capacity_ratio);

    let jet_thurst: f64 = (a * b * c).sqrt();
    let press_thrust: f64 = 
        (p_exhaust / p_chamber - p_atm / p_chamber) * expansion_ratio;
    let c_f: f64 = jet_thurst + press_thrust;
    return c_f
}


/// Calculate specific impulse of engine
/// 
/// Inputs
/// ------
/// heat_capacity_ratio : `f64`
///     Heat capacity ratio    
/// 
/// expansion_ratio : `f64`
///     Nozzle expansion ratio
/// 
/// t_chamber : `f64`
///     Chamber temperature
/// 
/// p_chamber : `f64`    
///     Chamber pressure
pub fn calc_engine_isp(
    heat_capacity_ratio: f64, 
    expansion_ratio: f64,
    t_chamber: f64,
    p_chamber: f64,
    p_atm: f64,
    p_exhaust: f64
) -> f64 {
    let thrust_coeff: f64 = calc_thrust_coeff(
        heat_capacity_ratio, 
        expansion_ratio, 
        p_chamber, 
        p_atm, 
        p_exhaust);
    let c_star: f64 = calc_characteristic_vel(heat_capacity_ratio, t_chamber);
    let isp: f64 = thrust_coeff * c_star;
    return isp
}
