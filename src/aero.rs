/*
Fluid dynamics
*/

#[path="./constants.rs"] mod constants;

pub fn calc_mass_flow(
    throat_area: f64,
    p_chamber: f64,
    c_star: f64
) -> f64 {
    let mass_flow: f64 = p_chamber * throat_area / c_star; 
    return mass_flow
}

pub fn calc_mach_number(
    velocity: f64,
    temperature: f64,
    heat_capacity_ratio: f64
) -> f64 {
    let speed_of_sound: f64 = 
        (heat_capacity_ratio * constants::GAS_CONST * temperature).sqrt();
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
    t_chamber: f64,
    molecular_weight: f64
) -> f64 {
    let a: f64 = 
        (2.0/(heat_capacity_ratio + 1.0)).powf(
            -(heat_capacity_ratio + 1.0)/(2.0*(heat_capacity_ratio - 1.0))); 
    let b: f64 = 
        ((constants::GAS_CONST * t_chamber) / 
        (heat_capacity_ratio * molecular_weight)).sqrt();
    let c_star: f64 = b * a;
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

pub fn calc_engine_isp(
    thrust_coeff: f64,
    c_star: f64
) -> f64 {
    let isp: f64 = thrust_coeff * c_star;
    return isp
}