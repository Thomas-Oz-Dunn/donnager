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
    r_1: f64,
    r_2: f64,
    v_0: f64
) -> f64 {
    let delta_v_1: f64 = v_0 * (((2.0 * r_2)/(r_1 + r_2)).sqrt()- 1.0);
    let delta_v_2: f64 = v_0 * ((r_1 / r_2).sqrt() + (2.0 * r_1)/(r_2 *(1.0 + r_2 / r_1)).sqrt());
    let delta_v_total: f64 = delta_v_1 + delta_v_2;
    delta_v_total
}

pub fn calc_stationary_orbit(
    mass_0: f64,
    period: f64
) -> f64 {
    let grav_param: f64 = mass_0 * constants::GRAV_CONST;
    let r: f64 = (grav_param.sqrt() * period / (2.0 * PI)).powf(2.0 / 3.0);
    r
}