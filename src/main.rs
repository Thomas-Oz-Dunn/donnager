/*
Orbital systems modelling in Rust

*/

use nalgebra as na;
use na::{Vector3};

use donnager as dgr;
use dgr::cosmos::orbit;
use dgr::propulsion::ballistics;
use dgr::constants;

fn main() {
    // Constants
    let earth: orbit::Body = orbit::Body {
        name: "Earth".to_string(),
        grav_param: constants::EARTH_GRAV_PARAM,
        eq_radius: constants::EARTH_RADIUS_EQUATOR,
        rotation_rate: constants::EARTH_ROT_RATE
    };
    
    // Inputs
    let n_stage: i32 = 1;
    let mass_1: f64 = 1.0; // kg
    let engine_isp: f64 = 300.0; // s
    let launch_pos_llh: Vector3<f64> = Vector3::new(28.396837, -80.605659, 0.0); // Cape Kennedy Lat Lon Height
    
    // Calculation
    let altitude: f64 = 408000.0;
    let radius: f64 = altitude + earth.eq_radius;
    let delta_v: f64 = earth.calc_orbital_velocity(radius);
    let surface_vel: f64 = earth.calc_surface_vel(
        launch_pos_llh);
    let net_delta_v: f64 = delta_v - surface_vel;
    let grav_acc: f64 = earth.grav_param / radius.powi(2);
    let mass_ratio: f64 = ballistics::calc_mass_ratio(net_delta_v, engine_isp, grav_acc);
    let mass_fuel: f64 = mass_1 * mass_ratio;    

    // Results
    println!("{} kg of fuel to get {} kg to {} m alt on {} stage", mass_fuel, mass_1, altitude, n_stage);

}
