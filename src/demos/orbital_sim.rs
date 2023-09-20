/*
Orbital systems modelling Application in Rust
*/
use std::{fs::File, io::{BufWriter, Write}};
use nalgebra::Vector3;

use donnager::donnager::{
    gravity as grav, 
    spacetime as xyzt, 
    constants as cst
};

fn main() {
    // TODO-TD: add cli parser

    let is_write_kml: bool = true;

    // Orbiting Body
    let tle_str: &str = "ISS
    1 25544U 98067A   23060.72453421  .00027779  00000-0  50068-3 0  9993
    2 25544  51.6419 140.8390 0005926  43.3718 100.9839 15.49547192385127";
    
    let orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_tle(tle_str.to_string());
    let frame: xyzt::ReferenceFrames = xyzt::ReferenceFrames::Perifocal;

    orbit.show(frame);

    let dt: f64 = 10.1;
    let new_orb: grav::kepler::Orbit = orbit.propogate(dt);
    let motion_ecef: Vector3<f64> = new_orb.calc_motion(
        0., 
        xyzt::ReferenceFrames::RotationalCartesian,
        0
    );

    // Earth
    let earth: xyzt::Body = xyzt::Body {
        name: String::from("Earth"),
        grav_param: cst::EARTH::GRAV_PARAM,
        eq_radius: cst::EARTH::RADIUS_EQUATOR,
        rotation_rate: cst::EARTH::ROT_RATE,
        sidereal_day_hours: cst::EARTH::SIDEREAL_DAY,
        eccentricity: cst::EARTH::SURFACE_ECC
    };

    let p_lla: Vector3<f64> = xyzt::ecef_to_lla(
        motion_ecef, 
        earth
    );
    
    // TODO: compare against observed baselines
    println!("{} is at {} altitude above {} deg N and {} deg E at {}",
        new_orb.name, p_lla[2], p_lla[0], p_lla[1], new_orb.epoch
    );

    // Plot over earth map

    // Write KML
    if is_write_kml{
        let mut kml_file = File::create("orbit.kml").unwrap();
        let mut kml_writer = BufWriter::new(&mut kml_file);
        let kml_str = orbit.to_kml();
        kml_writer.write_all(kml_str.as_bytes()).unwrap();
    }

}

