/*
Fuel calculations for launch
*/

use nalgebra::Vector3;
    
use donnager::donnager::{
    constants as cst, 
    spacetime as xyzt, 
    propulsion as prop, 
    dynamics as dynam
};


fn main() {
    
    // TODO-TD: add cli parsing

    // Earth
    let earth: xyzt::Body = xyzt::Body {
        name: String::from("Earth"),
        grav_param: cst::EARTH::GRAV_PARAM,
        eq_radius: cst::EARTH::RADIUS_EQUATOR,
        rotation_rate: cst::EARTH::ROT_RATE,
        sidereal_day_hours: cst::EARTH::SIDEREAL_DAY,
        eccentricity: cst::EARTH::SURFACE_ECC
    };

    // Config
    let launch_site: xyzt::SurfacePoint = xyzt::SurfacePoint {
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
        mass: 5.0,
        engine: payload_eng
    };

    let stage1_eng: prop::engine::Engine = prop::engine::Engine{
        name: "Hydrolox Dual Flow".to_string(),
        engine_type: prop::engine::EngineType::Chemical,
        isp: 300.0
    };

    let stage_1: dynam::vehicle::Vehicle = dynam::vehicle::Vehicle {
        name: "Stage_1".to_string(),
        mass: 250.0, // TODO: tune
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
    let mass_fuel: Vec<f64> = launch_vehicle.calc_mass_fuel(
        delta_v, 
        launch_site
    );

    // Results
    println!(
        "\n{:.4} kg of fuel to get {} kg to {} m alt", 
        mass_fuel[0], 
        launch_vehicle.stages[1].mass, 
        altitude
    );
    
}
