/*
Orbital systems modelling in Rust
*/

use nalgebra as na;
use na::Vector3;

use donnager::propulsion as prop;
use donnager::constants as cst;
use donnager::cosmos as cosm;
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

    let payload: dynam::vehicle::Vehicle = dynam::vehicle::Vehicle {
        name: "Satellite_1".to_string(),
        mass_0: 5.0,
        mass_prop: 0.001,
        engine_type: "Electric Ion".to_string(),
        engine_isp: 500.0,
    };

    let stage_1: dynam::vehicle::Vehicle = dynam::vehicle::Vehicle {
        name: "Stage_1".to_string(),
        mass_0: 250.0, // TODO: tune
        mass_prop: 500.0, // TODO: tune
        engine_type: "Hydrolox".to_string(),
        engine_isp: 300.0,
    };

    let launch_vehicle: dynam::vehicle::Multistage = dynam::vehicle::Multistage{
        name: "Launcher_7".to_string(),
        stages: [stage_1, payload].to_vec()
    };

    // Inputs
    let altitude: f64 = 408000.0;  // LEO

    // Calculation
    let delta_v: f64 = launch_site.calc_delta_v(altitude);
    let grav_acc: f64 = earth.calc_grav_acc(altitude + launch_site.calc_surface_radius());
    let mass_fuel: f64 = launch_vehicle.calc_mass_fuel(delta_v, grav_acc);
    
    // Results
    println!("\n{:.4} kg of fuel to get {} kg to {} m alt", mass_fuel, launch_vehicle.stages[1].mass_0 , altitude);
    
}
