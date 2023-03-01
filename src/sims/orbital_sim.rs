/*
Orbital systems modelling Application in Rust
*/

use nalgebra as na;
use na::Vector3;

use donnager::gravity as grav;

fn main() {

    let tle_str = 
    "UMBRA-2001
    1 48906U 21059AD  23060.23450352  .00060402  00000-0  27500-2 0  9997
    2 48906  97.5592 190.8803 0007011 262.3013  97.7426 15.20439300 92688";

    let orbit: grav::kepler::Orbit = Orbit.from_tle(tle_str);
    let eval_datetime: DateTime = DateTime();
    let p_eci: Vec<f64> = orbit.calc_position(eval_datetime);
    
    // TODO: compare against observed baselines

}