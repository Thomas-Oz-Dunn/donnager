/*
Calculate next overhead pass windows
*/


use nalgebra::Vector3;

use donnager::donnager::{
    gravity as grav,
    constants as cst, 
    spacetime as xyzt};


fn main() {
    // Orbiting Body
    let tle_str: &str = "ISS
    1 25544U 98067A   23060.72453421  .00027779  00000-0  50068-3 0  9993
    2 25544  51.6419 140.8390 0005926  43.3718 100.9839 15.49547192385127";
    
    let orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_tle(tle_str.to_string());
    let frame = xyzt::ReferenceFrames::Cartesian;

    let p_ecef = orbit.calc_motion(time_since_epoch, frame, 0);

    let ground_station: Vector3<f64> = vec![lat,   lon ,  height];

    


}