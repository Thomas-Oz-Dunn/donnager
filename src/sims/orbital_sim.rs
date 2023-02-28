/*
Orbital systems modelling Application in Rust
*/

use nalgebra as na;
use na::Vector3;

use donnager::gravity as grav;

fn main() {

    let tle_str = "";
    let orbit: grav::kepler::Orbit = Orbit.from_tle(tle_str);
    let eval_datetime: DateTime = DateTime();
    let p_eci: Vec<f64> = orbit.calc_position(eval_datetime);
    
    // TODO: compare against observed baselines

}