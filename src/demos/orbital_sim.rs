/*
Orbital systems modelling Application in Rust
*/

use nalgebra as na;
use na::Vector3;

use donnager::donnager::{gravity as grav, cosmos as cosm, constants as cst};

fn main() {

    let tle_str = "ISS
    1 25544U 98067A   23060.72453421  .00027779  00000-0  50068-3 0  9993
    2 25544  51.6419 140.8390 0005926  43.3718 100.9839 15.49547192385127";

    let orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_tle(tle_str.to_string());
    let frame = cosm::spacetime::ReferenceFrames::ECEF;
    let dt: f64 = 100.1;
    let new_orb = orbit.propogate(dt);
    let p_v_ecef: (Vector3<f64>, Vector3<f64>) = new_orb.calc_pos_vel(0., frame);

    let earth: cosm::spacetime::Body = cosm::spacetime::Body {
        name: "Earth".to_string(),
        grav_param: cst::EARTH_GRAV_PARAM,
        eq_radius: cst::EARTH_RADIUS_EQUATOR,
        rotation_rate: cst::EARTH_ROT_RATE,
        eccentricity: cst::EARTH_ECC
    };
    let p_lla: Vector3<f64> = earth.xyz_to_geodetic(p_v_ecef.0);
    
    // TODO: compare against observed baselines
    println!("{} is at {} altitude above {} deg N and {} deg E at {}",
        new_orb.name, p_lla[2], p_lla[0], p_lla[1], new_orb.epoch
    )
}
