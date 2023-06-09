/*
Interplanetary Mission Planner
*/
use std::f64::consts::PI;
use polars::prelude::*;
use chrono::{DateTime, Utc};
use nalgebra::{Vector3, Matrix3, self as na};

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
) -> Vector3<f64> {
    let frame = xyzt::ReferenceFrames::InertialCartesian;
    
    for launch_time in start_date_time.timestamp()..stop_date_time.timestamp() {
        let radius_1: Vector3<f64> = orbit_1.calc_motion(launch_time as f64, frame)[0];
        let radius_2: Vector3<f64> = orbit_2.calc_motion(launch_time as f64, frame)[0];
        let v_inf: f64 = calc_esc_vel(
            cst::SUN::GRAV_PARAM,
            radius_1.norm(), 
            radius_2.norm());
        }

    let lambert = lambert_solve(orb_dpt, orb_arr);
            
    let dv_dpt = man_lambert.impulses[0][1].norm();
    let dv_arr = man_lambert.impulses[1][1].norm();

    return [dv_dpt, dv_arrs]
}


/// Izzo based lambert solver
///
/// Inputs
/// ------
/// orbit_i: `kepler::Orbit`
///     Initial orbit
///
/// orbit_f: `kepler::Orbit`
///     Final orbit
///
/// Returns
/// -------
/// 

pub fn lambert_solve(
    orbit_i: kepler::Orbit,
    orbit_f: kepler::Orbit,
    prograde_sign: f64
){
    let k: f64 = orbit_i.central_body.grav_param;
    let r_i: Vector3<f64> = orbit_i.calc_motion(
        0.0, xyzt::ReferenceFrames::InertialCartesian)[0];
    let r_f: Vector3<f64> = orbit_f.calc_motion(
        0.0, xyzt::ReferenceFrames::InertialCartesian)[0];
    let tof = orbit_f.epoch - orbit_i.epoch;
    let chord: Vector3<f64> = r_f - r_i;

    let c_norm: f64 = chord.norm();
    let r_i_norm: f64 = r_i.norm();
    let r_f_norm: f64 = r_f.norm();
    let semi_perim: f64 = (r_i_norm + r_f_norm + c_norm) * 0.5;

    let u_r_i: Vector3<f64> = r_i / r_i_norm;
    let u_r_f: Vector3<f64> = r_f / r_f_norm;

    let i_h: Vector3<f64> = u_r_i.cross(&u_r_f);
    let u_h: Vector3<f64> = i_h / i_h.norm();

    let min: f64 = (1.0f64).min(c_norm / semi_perim);
    let mut ll: f64 = prograde_sign * (1.0 - min).powf(0.5);
    if u_h[2] < 0.0 {
        ll: f64 = -ll;
        let u_t_i: Vector3<f64> = prograde_sign * u_r_i.cross(&u_h);
        let u_t_f: Vector3<f64> = prograde_sign * u_r_f.cross(&u_h);
    }
    else{
        let u_t_i: Vector3<f64> = prograde_sign * u_h.cross(&u_r_i);
        let u_t_f: Vector3<f64> = prograde_sign * u_h.cross(&u_r_f);
    }

    let time = tof * ((2. * k / semi_perim.powi(3)).powf(0.5) as i32);

    
    let M_max: f64 = (time / PI).floor();
    let T_00: f64 = ll.acos() + ll * (1. - ll.powi(2)).powf(0.5);  // T_xM
    
    // FIXME-TD: Translate into Rust -V
    // # Refine maximum number of revolutions if necessary
    // if T < T_00 + M_max * pi and M_max > 0:
    //     _, T_min = _compute_T_min(ll, M_max, numiter, rtol)
    //     if T < T_min:
    //         M_max -= 1

    // # Check if a feasible solution exist for the given number of revolutions
    // # This departs from the original paper in that we do not compute all solutions
    // if M > M_max:
    //     raise ValueError("No feasible solution, try lower M")

    // # Initial guess
    // x_0 = _initial_guess(T, ll, M, lowpath)

    // # Start Householder iterations from x_0 and find x, y
    // x = _householder(x_0, T, ll, M, rtol, numiter)
    // for ii in range(maxiter):
    //     y = np.sqrt(1 - ll**2 * (1 - x_0**2)) 
    //     if M == 0 and np.sqrt(0.6) < x < np.sqrt(1.4):
    //         eta = y - ll * x_0
    //         S_1 = (1 - ll - x_0 * eta) * 0.5
    //         Q = 4 / 3 * hyp2f1b(S_1)
    //         T_ = (eta**3 * Q + 4 * ll * eta) * 0.5
    //     else:
    //         if -1 <= x < 1:
    //             psi = np.arccos(x * y + ll * (1 - x**2))
    //         elif x > 1:
    //             psi = np.arcsinh((y - x * ll) * np.sqrt(x**2 - 1))
    //         else:
    //             psi = 0.0

    //         psi = _compute_psi(x_0, y, ll)
    //         T_ = np.divide(
    //             np.divide(psi + M * pi, np.sqrt(np.abs(1 - x_0**2))) - x_0 + ll * y,
    //             (1 - x_0**2),
    //         )

    //     fval =  T_ - T0
        
    //     T = fval + T0
    //     fder = (3 * T * x_0 - 2 + 2 * ll**3 * x_0 / y) / (1 - x_0**2) 
    //     fder2 = (3 * T + 5 * x_0 * dT + 2 * (1 - ll**2) * ll**3 / y**3) / (1 - x**2)
    //     fder3 = (
    //         7 * x_0 * ddT + 8 * dT - 6 * (1 - ll**2) * ll**5 * x / y**5
    //     ) / (1 - x_0**2)

    //     # Householder step (quartic)
    //     p = x_0 - fval * (
    //         (fder**2 - fval * fder2 / 2) / (fder * (fder**2 - fval * fder2) + fder3 * fval**2 / 6)
    //     )

    //     if abs(x - x_0) < tol:
    //         return x
    //     x = p

    let y: f64 = (1. - ll.powi(2) * (1. - x.powi(2))).powf(0.5);

    let gamma: f64 = (k * semi_perim / 2.).powf(0.5);
    let rho: f64 = (r_i_norm - r_f_norm) / c_norm;
    let sigma: f64 = (1. - rho.powi(2)).powf(0.5);

    let v_r_i = gamma * ((ll * y - x) - rho * (ll * y + x)) / r_i;
    let v_r_f = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r_f;
    let v_t_i = gamma * sigma * (y + ll * x) / r_i;
    let v_t_f = gamma * sigma * (y + ll * x) / r_f;

    let v_i: Vector3<f64> = v_r_i * (r_i / r_i_norm) + v_t_i * i_t_i;
    let v_f: Vector3<f64> = v_r_f * (r_f / r_f_norm) + v_t_f * i_t_f;

    return (v_i, v_f)
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
