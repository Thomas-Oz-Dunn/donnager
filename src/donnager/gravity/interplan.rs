/*
Interplanetary Planner
*/


use crate::donnager::{spacetime as xyzt, constants as cst, gravity::kepler as kepler};


/// Calculate interplanetary mission escape velocity
/// 
/// Inputs
/// ------
/// orb_radius_0: `f64`
///     Starting planet heliocentric radius
/// 
/// orb_radius_f: `f64`
///     Ending planet heliocentric radius
pub fn calc_esc_vel(
    orb_radius_0: f64,
    orb_radius_f: f64,
) -> f64 {
    return (2. * cst::SUN::grav_param / orb_radius_0 - 2. * cst::SUN::grav_param /(orb_radius_0+orb_radius_f)).sqrt();
}


/// Calculate next hohmann transfer launch windows
/// 
pub fn calc_next_hohmann_launch_windows(
    start_datetime: DateTime<UTC>,
    orbit_1: kepler::Orbit,
    orbit_2: kepler::Orbit,
    n_windows: usize
) -> Vec<DateTime<UTC>>{

    let period_1 = orbit_1.calc_period();
    let period_2 = orbit_2.calc_period();
    let synodic_period = period_1 / period_2;

    

}




