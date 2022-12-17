/*
Orbital systems modelling in Rust

*/
use nalgebra as na;
use na::{Vector3};

mod donnager;
use donnager::{constants, ballistics, astro};

fn main() {
    let equatorial_radius: f64 = constants::EARTH_RADIUS_EQUATOR; 
    let grav_param: f64 = constants::EARTH_GRAV_PARAM;
    let rotation_rate: f64 = constants::EARTH_ROT_RATE;

    let mass_1: f64 = 1.0; // kg
    let altitude: f64 = 3.0e6; // 300km, LEO
    let radius_f: f64 = equatorial_radius + altitude; // m
    let engine_isp: f64 = 300.0; // s
    let pos_llh: Vector3<f64> = Vector3::new(28.396837, -80.605659, 0.0); // Cape Kennedy Lat Lon Height

    let delta_v: f64 = astro::calc_orbital_velocity(grav_param, radius_f);
    let surface_vel: f64 = astro::calc_surface_vel(
        rotation_rate, 
        equatorial_radius, 
        pos_llh);
    let net_delta_v: f64 = delta_v - surface_vel;
    let grav_acc: f64 = grav_param / radius_f.powi(2);
    let mass_ratio: f64 = ballistics::calc_mass_ratio(net_delta_v, engine_isp, grav_acc);
    let mass_fuel: f64 = mass_1 * mass_ratio;    
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