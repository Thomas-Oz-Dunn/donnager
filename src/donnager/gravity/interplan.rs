/*
Interplanetary Mission Planner
*/
use std::f64::consts::PI;
use na::Vector;
use chrono::{DateTime, Utc};
use nalgebra::{Vector3, Matrix3, self as na};

use crate::donnager::{
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
) -> Vec<f64> {
    let frame = xyzt::ReferenceFrames::InertialCartesian;
    
    for launch_time in start_date_time.timestamp()..stop_date_time.timestamp() {
        let radius_1: Vector3<f64> = orbit_1.calc_motion(launch_time as f64, frame)[0];
        let radius_2: Vector3<f64> = orbit_2.calc_motion(launch_time as f64, frame)[0];

        let lambert_vels = lambert_solve(
            orbit_1, 
            orbit_2, 
            1.0);
                
        let dv_dpt = lambert_vels[0] - orbit_1.calc_motion(start_date_time as f64, frame)[1];
        let dv_arr = lambert_vels[1] - orbit_1.calc_motion(stop_date_time as f64, frame)[1];
        let dv_tot = dv_dpt.norm() + dpt_arr.norm();
        }
    return dv_tot
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
/// prograde_sign: `f64`
///     Prograde vs retrograde
/// 
/// Returns
/// -------
/// 
pub fn lambert_solve(
    orbit_i: kepler::Orbit,
    orbit_f: kepler::Orbit,
    prograde_sign: f64
) -> (Vector3<f64>, Vector3<f64>) {
    let k: f64 = orbit_i.central_body.grav_param;
    let r_i: Vector3<f64> = orbit_i.calc_motion(
        0.0, xyzt::ReferenceFrames::InertialCartesian)[0];
    let r_f: Vector3<f64> = orbit_f.calc_motion(
        0.0, xyzt::ReferenceFrames::InertialCartesian)[0];
    let tof = orbit_f.epoch - orbit_i.epoch;
    let chord: Vector3<f64> = r_f - r_i;

    // Normalize operations
    let c_norm: f64 = chord.norm();
    let r_i_norm: f64 = r_i.norm();
    let r_f_norm: f64 = r_f.norm();
    let semi_perim: f64 = (r_i_norm + r_f_norm + c_norm) * 0.5;

    let u_r_i: Vector3<f64> = r_i / r_i_norm;
    let u_r_f: Vector3<f64> = r_f / r_f_norm;

    let i_h: Vector3<f64> = u_r_i.cross(&u_r_f);
    let u_h: Vector3<f64> = i_h / i_h.norm();

    let min: f64 = (1.0f64).min(c_norm / semi_perim);
    let mut ll: f64 = prograde_sign * (1.0 - min).sqrt();
    if u_h[2] < 0.0 {
        ll = -ll;
        let u_t_i: Vector3<f64> = prograde_sign * u_r_i.cross(&u_h);
        let u_t_f: Vector3<f64> = prograde_sign * u_r_f.cross(&u_h);
    }
    else{
        let u_t_i: Vector3<f64> = prograde_sign * u_h.cross(&u_r_i);
        let u_t_f: Vector3<f64> = prograde_sign * u_h.cross(&u_r_f);
    }

    // Dimensionless time parameters
    let time = tof * ((2. * k / semi_perim.powi(3)).sqrt() as i32);
    
    let M_max: f64 = time / (PI as i32);
    let T_00: f64 = ll.acos() + ll * (1. - ll.powi(2)).sqrt();  // T_xM
    
    // FIXME-TD: Translate into Rust -V
    // # Refine maximum number of revolutions if necessary
    if (time < T_00 + M_max * PI) && (M_max > 0){
        //     _, T_min = _compute_T_min(ll, M_max, numiter, rtol)
        //     if T < T_min:
        //         M_max -= 1

    }

    // # Initial guess
    // x_0 = _initial_guess(T, ll, M, lowpath)

    // # Start Householder iterations from x_0 and find x, y
    // x = _householder(x_0, T, ll, M, rtol, numiter)

    // for ii in 0..maxiter {
    //     let y = (1 - ll.powi(2) * (1 - x_0.powi(2))).sqrt();

    //     if M == 0.0 && (0.6).sqrt() < x && x < (1.4).sqrt(){
    //         let mut eta = y - ll * x_0;
    //         S_1 = (1 - ll - x_0 * eta) * 0.5
    //         Q = 4 / 3 * hyp2f1b(S_1)
    //         T_= (eta**3 * Q + 4 * ll * eta) * 0.5
    //     } else {

    //         if -1 <= x < 1{
    //             psi = np.arccos(x * y + ll * (1 - x**2));
    //         }
    //         else if (x > 1){
    //             psi = np.arcsinh((y - x * ll) * (x**2 - 1).sqrt());
    //         }
    //         else{
    //             psi = 0.0;
    //         }

    //         psi = _compute_psi(x_0, y, ll)
    //         T_ = np.divide(
    //             np.divide(psi + M * pi, np.sqrt(np.abs(1 - x_0**2))) - x_0 + ll * y,
    //             (1 - x_0**2),
    //         )
    //     }

    //     fval =  T_ - T0
        
    //     T = fval + T0
    //     fder = (3 * T * x_0 - 2 + 2 * ll**3 * x_0 / y) / (1 - x_0**2) 
    //     fder2 = (3 * T + 5 * x_0 * dT + 2 * (1 - ll**2) * ll**3 / y**3) / (1 - x**2)
    //     fder3 = (
    //         7 * x_0 * ddT + 8 * dT - 6 * (1 - ll**2) * ll**5 * x / y**5
    //     ) / (1 - x_0**2)

    //     // # Householder step (quartic)
    //     p = x_0 - fval * (
    //         (fder**2 - fval * fder2 / 2) / (fder * (fder**2 - fval * fder2) + fder3 * fval**2 / 6)
    //     )

    //     if abs(x - x_0) < tol:
    //         return x
    //     x = p
    // }

    let y: f64 = (1. - ll.powi(2) * (1. - x.powi(2))).sqrt();

    let gamma: f64 = (k * semi_perim / 2.).sqrt();
    let rho: f64 = (r_i_norm - r_f_norm) / c_norm;
    let sigma: f64 = (1. - rho.powi(2)).sqrt();

    // Radial and Tangential Components
    let v_r_i = gamma * ((ll * y - x) - rho * (ll * y + x)) / r_i;
    let v_r_f = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r_f;
    let v_t_i = gamma * sigma * (y + ll * x) / r_i;
    let v_t_f = gamma * sigma * (y + ll * x) / r_f;

    // Velocity vectors
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
