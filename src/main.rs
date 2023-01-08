/*
Orbital systems modelling in Rust
*/

use nalgebra as na;
use na::Vector3;

use donnager::propulsion as prop;
use donnager::constants as cst;
use donnager::cosmos as cosm;

fn main() {
    // Config
    let earth: cosm::objects::Body = cosm::objects::Body {
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

    let payload: cosm::objects::Vehicle = cosm::objects::Vehicle {
        name: "Satellite_1".to_string(),
        mass_0: 5.0,
        mass_prop: 0.001,
        engine_type: "Electric Ion".to_string(),
        engine_isp: 500.0,
    };

    // Inputs
    let n_stage: i32 = 1;
    let launch_engine_isp: f64 = 300.0; // s
    let altitude: f64 = 408000.0;  // LEO

    // Calculation
    let delta_v: f64 = launch_site.calc_delta_v(altitude);
    let grav_acc: f64 = earth.calc_grav_acc(altitude + launch_site.calc_surface_radius());
    let mass_ratio: f64 = prop::ballistics::calc_mass_ratio(delta_v, launch_engine_isp, grav_acc);
    let mass_fuel: f64 = (payload.mass_0 + payload.mass_prop) * mass_ratio;    

    // Results
    println!("\n{:.4} kg of fuel to get {} kg to {} m alt on {} stage", mass_fuel, payload.mass_0 , altitude, n_stage);

}
