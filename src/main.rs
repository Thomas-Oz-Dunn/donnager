/*
Orbital systems modelling in Rust
*/

use nalgebra as na;
use na::Vector3;

use donnager::constants as cst;
use donnager::cosmos as cosm;
use donnager::propulsion as prop;
use donnager::dynamics as dynam;


fn main() {
    // Config
    let earth: cosm::grav::Body = cosm::grav::Body {
        name: "Earth".to_string(),
        grav_param: cst::EARTH_GRAV_PARAM,
        eq_radius: cst::EARTH_RADIUS_EQUATOR,
        rotation_rate: cst::EARTH_ROT_RATE,
        eccentricity: cst::EARTH_ECC
    };
    
    let launch_site: cosm::space::SurfacePoint = cosm::space::SurfacePoint {
        name: "Cape Canaveral Launch Site".to_string(),
        body: earth.clone(),
        pos_lla: Vector3::new(28.396837, -80.605659, 0.0)
    };

    let payload_eng: prop::engine::Engine = prop::engine::Engine{
        name: "Payload_Engine".to_string(),
        engine_type: prop::engine::EngineType::Electric,
        isp: 500.0
    };

    let payload: dynam::vehicle::Vehicle = dynam::vehicle::Vehicle {
        name: "Satellite_1".to_string(),
        mass_0: 5.0,
        engine: payload_eng
    };

    let stage1_eng: prop::engine::Engine = prop::engine::Engine{
        name: "Hydrolox Dual Flow".to_string(),
        engine_type: prop::engine::EngineType::Chemical,
        isp: 300.0
    };

    let stage_1: dynam::vehicle::Vehicle = dynam::vehicle::Vehicle {
        name: "Stage_1".to_string(),
        mass_0: 250.0, // TODO: tune
        engine: stage1_eng
    };

    let launch_vehicle: dynam::vehicle::Multistage = dynam::vehicle::Multistage{
        name: "Launcher_7".to_string(),
        stages: [stage_1, payload].to_vec()
    };

    // Inputs
    let altitude: f64 = 408000.0;  // LEO

    // Calculation
    let delta_v: f64 = launch_site.calc_delta_v(altitude);
    let mass_fuel: Vec<f64> = launch_vehicle.calc_mass_fuel(delta_v, launch_site);

    // Results
    println!("\n{:.4} kg of fuel to get {} kg to {} m alt", mass_fuel[0], launch_vehicle.stages[1].mass_0 , altitude);
    
}
