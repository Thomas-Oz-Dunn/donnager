/*
Orbital systems modelling in Rust
*/

use nalgebra as na;
use na::{Vector3};

use donnager as dgr;
use dgr::*;

fn main() {
    // Constants
    let earth: cosmos::objects::Body = cosmos::objects::Body {
        name: "Earth".to_string(),
        grav_param: constants::EARTH_GRAV_PARAM,
        eq_radius: constants::EARTH_RADIUS_EQUATOR,
        rotation_rate: constants::EARTH_ROT_RATE,
        eccentricity: constants::EARTH_ECC
    };
    
    let launch_site: cosmos::space::SurfacePoint = cosmos::space::SurfacePoint {
        name: "Cape Canaveral Launch Site".to_string(),
        body: earth.clone(),
        pos_lla: Vector3::new(28.396837, -80.605659, 0.0)
    };

    let sat_1: cosmos::objects::Craft = cosmos::objects::Craft {
        name: "Satellite_1".to_string(),
        mass: 5.0,
        engine_type: "Electric Ion".to_string(),
        engine_isp: 500.0,
    };

    // Inputs
    let n_stage: i32 = 1;
    let launch_engine_isp: f64 = 300.0; // s
    let altitude: f64 = 408000.0;

    // Calculation
    let radius: f64 = altitude + earth.eq_radius;
    let delta_v: f64 = earth.calc_orbital_velocity(radius);
    let surface_vel: f64 = launch_site.calc_surface_vel();
    let net_delta_v: f64 = delta_v - surface_vel;
    let grav_acc: f64 = earth.calc_grav_acc(radius);
    let mass_ratio: f64 = propulsion::ballistics::calc_mass_ratio(net_delta_v, launch_engine_isp, grav_acc);
    let mass_fuel: f64 = sat_1.mass * mass_ratio;    

    // Results
    println!("\n{:.4} kg of fuel to get {} kg to {} m alt on {} stage", mass_fuel, sat_1.mass , altitude, n_stage);

}
