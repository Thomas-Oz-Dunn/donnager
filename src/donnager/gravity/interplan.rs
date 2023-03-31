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

    let epoch_time = start_datetime.timestamp as f64 - orbit_1.epoch.timestamp as f64;;
    let true_anonmaly_0_1 = orbit_1.calc_true_anomaly(epoch_time);

    let epoch_time = start_datetime.timestamp as f64 - orbit_2.epoch.timestamp as f64;;
    let true_anonmaly_0_2 = orbit_2.calc_true_anomaly(epoch_time);

}



/// Calculate lagrange point locations in DU
/// 
/// Inputs
/// ------
/// mass_1: `f64`
///     Larger mass
/// 
/// mass_2: `f64`
///     Smaller mass
pub fn calc_lagrange_points(
    mass_1: f64,
    mass_2: f64
) -> Vec<Vec3<f64>> {
    let mass_ratio: f64= mass_2 / (mass_1 + mass2);
    let l1 = Vec3::new([xl1,0.,0.]);
    let l2 = Vec3::new([xl2 ,0.,0.]);
    let l3 = Vec3::new([xl3 ,0.,0.]);
    let l4 = Vec3::new([mass_ratio - 0.5, -(3.).sqrt()/2., 0.]);
    let l5 = Vec3::new([mass_ratio - 0.5, (3.).sqrt()/2., 0.]);
    return vec![l1, l2, l3, l4, l5]
}


#[cfg(test)]
mod interplan_tests {

}
