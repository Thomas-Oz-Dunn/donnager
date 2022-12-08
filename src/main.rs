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
    let engine_isp: f64 = 300.0; // s

    let delta_v: f64 = astro::calc_orbital_velocity(mass_0, radius_f);
    let grav_acc: f64 = astro::calc_grav_acc(mass_0, radius_0);
    let mass_ratio: f64 = ballistics::calc_mass_ratio(delta_v, engine_isp, grav_acc);
    let mass_fuel: f64 = mass_1 * mass_ratio;    
    println!("{} kg of fuel", mass_fuel);
}

/*
TODO-TD: 
Convert scalar pos, vel, acc into vectors.
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