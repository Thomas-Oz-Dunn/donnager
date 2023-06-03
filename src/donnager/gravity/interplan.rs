/*
Interplanetary Mission Planner
*/
use polars::prelude::*;
use chrono::{DateTime, Utc};
use plotters::prelude::*;
use std::ops::Range;

use crate::donnager::{
    constants as cst, 
    gravity::kepler as kepler, 
    spacetime as xyzt};

/// Calculate interplanestary mission escape velocity
/// 
/// Inputs
/// ------
/// orb_radius_0: `f64`
///     Starting planet heliocentric radius
/// 
/// orb_radius_f: `f64`
///     Ending planet heliocentric radius
pub fn calc_esc_vel(
    grav_param: f64,
    orb_radius_0: f64,
    orb_radius_f: f64,
) -> f64 {
    let mean_radius: f64 = (orb_radius_0 + orb_radius_f) / 2.;
    return (grav_param * (2. / orb_radius_0 - 1. / mean_radius)).sqrt();
}

/// Get ephemeris dataframe
pub fn get_ephemeris() -> DataFrame {

    let ephemeris = df! (
        "number" => &[1, 2, 3, 4, 5, 6, 7, 8],
        "names" => &["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"],
        "groups" => &["Inner", "Inner", "Inner", "Inner", "Outer", "Outer", "Outer", "Outer"],
        "type" => &["Rocky", "Rocky", "Rocky", "Rocky", "Gas Giant", "Gas Giant", "Ice Giant", "Ice Giant"],
        "mass" => &[
            cst::MERCURY::MASS, 
            cst::VENUS::MASS, 
            cst::EARTH::MASS, 
            cst::MARS::MASS, 
            cst::JUPITER::MASS, 
            cst::SATURN::MASS, 
            cst::URANUS::MASS, 
            cst::NEPTUNE::MASS],

        "orbit_semi_major" => &[
            cst::MercurySunOrbit::SEMI_MAJOR, 
            cst::VenusSunOrbit::SEMI_MAJOR, 
            cst::EarthSunOrbit::SEMI_MAJOR, 
            cst::MarsSunOrbit::SEMI_MAJOR, 
            None, 
            None, 
            None, 
            None],

        "orbit_ecc" => &[
            cst::MercurySunOrbit::ECC, 
            cst::VenusSunOrbit::ECC, 
            cst::EarthSunOrbit::ECC, 
            cst::MarsSunOrbit::ECC, 
            None, 
            None, 
            None, 
            None],
        
        "argument_perigelion" => &[
            cst::MercurySunOrbit::ARG_PERIHELION, 
            cst::VenusSunOrbit::ARG_PERIHELION, 
            cst::EarthSunOrbit::ARG_PERIHELION, 
            cst::MarsSunOrbit::ARG_PERIHELION, 
            None, 
            None, 
            None, 
            None],

    ).unwrap();

    return ephemeris
}


/// Calculate orbital transfer porkchop plots
/// 
/// Inputs
/// ------
/// datetime_launch_window: `(DateTime<Utc>, DateTime<Utc>)`
///     Start and end datetime for plot
/// 
/// orbit_1: `Orbit`
///     Orbit of starting planet
/// 
/// orbit_2: `Orbit`
///     Orbit of ending planet
pub fn calc_porkchop_plots(
    start_date_time: DateTime<Utc>,
    stop_date_time: DateTime<Utc>,
    orbit_1: kepler::Orbit,
    orbit_2: kepler::Orbit
) -> Vector<f64> {
    let frame = xyzt::ReferenceFrames::InertialCartesian;
    
    for launch_time in start_date_time.timestamp()..stop_date_time.timestamp() {
        let motion1 = orbit_1.calc_motion(launch_time as f64, frame);
        let motion2 = orbit_2.calc_motion(launch_time as f64, frame);
        let v_inf: f64 = calc_esc_vel(
            cst::SUN::GRAV_PARAM,
            motion1[0].norm(), 
            motion2[0].norm());
        }

    let lambert = lambert(orb_dpt, orb_arr);
            
    // Get norm delta velocities
    let dv_dpt = nalgebra.norm(man_lambert.impulses[0][1]);
    let dv_arr = nalgebra.norm(man_lambert.impulses[1][1]);
            
    // Compute all the output variables
    let c3_launch = dv_dpt**2;
    let c3_arrival = dv_arr**2;

    return [dv_dpt, dv_arr, c3_launch, c3_arrival]
}


pub fn lambert_solve(
    orbit_1: kepler::Orbit,
    orbit_2: kepler::Orbit

) {
    let k = orbit_i.Body.GRAV_PARAM;
    let r_i = orbit_i.radial_distance;
    let r_f = orbit_f.r;
    let tof = orbit_f.epoch - orbit_i.epoch;
    let chord = r_f - r_i;

    let c_norm = nalgebra::norm(chord);
    let r_i_norm = nalgebra::norm(r_i);
    let r_f_norm = nalgebra::norm(r_f);
    let semi_perim = (r1_norm + r2_norm + c_norm) * 0.5;

    let i_r1 = r_i / r_i_norm;
    let i_r2 = r_f / r_f_norm;

    let i_h = nalgebra::cross(i_r1, i_r2);
    let i_h = i_h / nalgebra::norm(i_h);
    
    if i_h[2] < 0.0{
        let ll = -nalgebra.sqrt(1 - min(1.0, c_norm / semi_perim));
        let i_t_i = cross(i_r_i, i_h);
        let i_t_f = cross(i_r_f, i_h);
    }
    else{
        let ll = nalgebra.sqrt(1 - min(1.0, c_norm / semi_perim));
        let i_t_i = cross(i_h, i_r_i);
        let i_t_f = cross(i_h, i_r_f);
    }

    ll, i_t1, i_t2 = (ll, i_t1, i_t2) if prograde else (-ll, -i_t1, -i_t2)

    T = np.sqrt(2 * k / s**3) * tof

    x, y = _find_xy(ll, T, M, numiter, lowpath, rtol)

    gamma = np.sqrt(k * s / 2)
    rho = (r1_norm - r2_norm) / c_norm
    sigma = np.sqrt(1 - rho**2)

    V_r1, V_r2, V_t1, V_t2 = _reconstruct(
        x, y, r1_norm, r2_norm, ll, gamma, rho, sigma
    )

    v1 = V_r1 * (r1 / r1_norm) + V_t1 * i_t1
    v2 = V_r2 * (r2 / r2_norm) + V_t2 * i_t2

    return 
}

#[cfg(test)]
mod interplan_tests {

    use super::*;

    #[test]
    fn test_calc_esc_vel()
    {
        let grav_param: f64 = cst::SUN::GRAV_PARAM;
        let orb_radius_0: f64 = 1000000.;
        let orb_radius_f: f64 = 2000000.;
        let esc_vel: f64 = calc_esc_vel(grav_param, orb_radius_0, orb_radius_f);
        assert_eq!(esc_vel, 13302555.410020536);
    }


}
