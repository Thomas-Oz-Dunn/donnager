/*
Orbital systems modelling in Rust

*/
mod aero;
mod astro;
mod constants;
mod ballistics;

fn main() {
    let mass_0: f64 = constants::EARTH_MASS;
    let radius_0: f64 = constants::EARTH_RADIUS_EQUATOR; 

    let mass_1: f64 = 1.0; // kg
    let altitude: f64 = 3.0e6; // 300km, LEO
    let radius_f: f64 = radius_0 + altitude; // m
    let mut engine_isp: f64 = 300.0; // ss

    let c_star: f64 = aero::calc_characteristic_vel(
        heat_capacity_ratio, 
        t_chamber, 
        molecular_weight)
    let thrust_coeff: f64 = aero::calc_thrust_coeff(
        heat_capacity_ratio, 
        expansion_ratio, 
        p_chamber, 
        p_atm, 
        p_exhaust)
    engine_isp = aero::calc_engine_isp(thrust_coeff, c_star);

    let grav_param: f64 = mass_0 * constants::GRAV_CONST;
    let delta_v: f64 = astro::calc_orbital_velocity(grav_param, radius_f);
    let grav_acc: f64 = grav_param / radius_f.powi(2);
    let mass_ratio: f64 = ballistics::calc_mass_ratio(delta_v, engine_isp, grav_acc);
    let mass_fuel: f64 = mass_1 * mass_ratio;    
    println!("{} kg of fuel", mass_fuel);
}

/*
TODO-TD: 
N-stage trade study plots
TLE ingest
Trans Lunar Injections
Multithreading, Cloud Compute?
RK45 propogator
J2 perturbation
Launch cost calculator
Solar System Mineralogical data base query
Interplanetary Mission Plan (optimal launch windows, porkchop plot)
$ / kg (mineral X)
*/