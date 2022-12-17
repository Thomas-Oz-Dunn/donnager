/*
Orbital systems modelling in Rust

*/
use nalgebra as na;
use na::{Vector3};
use std::f64::consts::PI;

mod donnager;
use donnager::{constants, ballistics, astro};

fn main() {
    // Constants
    let equatorial_radius: f64 = constants::EARTH_RADIUS_EQUATOR; 
    let grav_param: f64 = constants::EARTH_GRAV_PARAM;
    let rotation_rate: f64 = constants::EARTH_ROT_RATE; // rad / s
    
    // Inputs
    let mass_1: f64 = 1.0; // kg
    let engine_isp: f64 = 300.0; // s
    let launch_pos_llh: Vector3<f64> = Vector3::new(28.396837, -80.605659, 0.0); // Cape Kennedy Lat Lon Height
    
    // Calculation
    let period: f64 = 2.0 * PI / rotation_rate; // s / 2* pi rad 
    let gso_radius: f64 = astro::calc_stationary_orbit(grav_param, period);
    let altitude: f64 = gso_radius - equatorial_radius;
    let delta_v: f64 = astro::calc_orbital_velocity(grav_param, gso_radius);
    let surface_vel: f64 = astro::calc_surface_vel(
        rotation_rate, 
        equatorial_radius, 
        launch_pos_llh);
    let net_delta_v: f64 = delta_v - surface_vel;
    let grav_acc: f64 = grav_param / gso_radius.powi(2);
    let mass_ratio: f64 = ballistics::calc_mass_ratio(net_delta_v, engine_isp, grav_acc);
    let mass_fuel: f64 = mass_1 * mass_ratio;    

    // Results
    println!("{} kg of fuel to get {} kg to {} m alt", mass_fuel, mass_1, altitude);

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
Comparitive propulsion techniques
    Liquid
        Monoprop
        Biprop
    Solid
    Electric
    Nuclear Thermal
Solar System Mineralogical data base query
Interplanetary Mission Plan (optimal launch windows, porkchop plot)
$ / kg (mineral X)
*/