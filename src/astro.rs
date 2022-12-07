use std::f64::consts::PI;
#[path="./constants.rs"] mod constants;

pub fn calc_grav_acc(
    mass_0: f64,
    radius: f64
) -> f64{
    // Calculate gravitational acceleration
    // 
    // Parameters
    // ---------
    // mass_0: f64
    //     Mass of object
    // 
    // radius: f64
    //     Radial distance to object
    // 
    // Returns
    // -------
    // grav_acc
    let grav_acc: f64 = constants::GRAV_CONST * mass_0 / radius.powi(2);
    grav_acc
}

pub fn calc_orbital_velocity(
    mass_0: f64,
    radius: f64
) -> f64{
    let energy: f64 = 2.0 * constants::GRAV_CONST * mass_0 / radius;
    let velocity: f64 = energy.sqrt();
    velocity
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
    delta_v_total
}

pub fn calc_stationary_orbit(
    mass: f64,
    period: f64
) -> f64 {
    let grav_param: f64 = mass * constants::GRAV_CONST;
    let r: f64 = (grav_param.sqrt() * period / (2.0 * PI)).powf(2.0 / 3.0);
    r
}