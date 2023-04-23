/*
Orbital systems modelling Application in Rust
*/
use std::{fs::File, io::{BufWriter, Write}};
use nalgebra::Vector3;

use donnager::donnager::{gravity as grav, spacetime as xyzt};

fn main() {

    let is_write_kml: bool = true;

    // Orbiting Body
    let tle_str: &str = "ISS
    1 25544U 98067A   23060.72453421  .00027779  00000-0  50068-3 0  9993
    2 25544  51.6419 140.8390 0005926  43.3718 100.9839 15.49547192385127";
    
    let orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_tle(tle_str.to_string());
    let frame = xyzt::ReferenceFrames::PFCL;

    orbit.show(frame);

    let dt: f64 = 10.1;
    let new_orb: grav::kepler::Orbit = orbit.propogate(dt);
    let p_v_ecef = new_orb.calc_pos_vel(0., xyzt::ReferenceFrames::ECEF);

    let p_lla: Vector3<f64> = xyzt::ecef_to_lla(p_v_ecef[0]);
    
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

