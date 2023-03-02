/*
Orbital systems modelling Application in Rust
*/

// use nalgebra as na;
// use na::Vector3;s

use donnager::gravity as grav;

fn main() {

    let tle_str = "ISS
    1 25544U 98067A   23060.72453421  .00027779  00000-0  50068-3 0  9993
    2 25544  51.6419 140.8390 0005926  43.3718 100.9839 15.49547192385127";

    let orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_tle(tle_str.to_string());
    println!("{:?}", orbit);
    // let eval_datetime: DateTime = DateTime();
    // let p_eci: Vec<f64> = orbit.calc_position(eval_datetime);
    
    // TODO: compare against observed baselines

}