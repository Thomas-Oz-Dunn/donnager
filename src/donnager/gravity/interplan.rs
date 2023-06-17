/*
Interplanetary Mission Planner
*/
use std::f64::consts::PI;
use chrono::{DateTime, Utc, Timelike};
use nalgebra::{Vector3};

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

/// Calculate delta v for mission
/// 
/// Inputs
/// ------
/// start_date_time: `DateTime<Utc>`
///     Start datetime for mission
/// 
/// stop_date_time: `DateTime<Utc>`
///     stop datetime for mission
/// 
/// orbit_1: `Orbit`
///     Orbit of starting planet
/// 
/// orbit_2: `Orbit`
///     Orbit of ending planet
pub fn calc_mission_delta_v(
    start_date_time: DateTime<Utc>,
    stop_date_time: DateTime<Utc>,
    orbit_i: kepler::Orbit,
    orbit_f: kepler::Orbit
) -> f64 {
    // TODO-TD: parallelize across start stop datetimes using rayon
    let grav_param: f64 = orbit_i.central_body.grav_param;
    let r_i: Vector3<f64> = orbit_i.calc_motion(
        0.0, xyzt::ReferenceFrames::InertialCartesian)[0];
    let r_f: Vector3<f64> = orbit_f.calc_motion(
        0.0, xyzt::ReferenceFrames::InertialCartesian)[0];
    let tof = (orbit_f.epoch - orbit_i.epoch).num_seconds() as f64;
    let (dv1, dv2) = lambert_solve(
        grav_param,  
        r_i,  
        r_f, 
        tof, 
        1.0);
        
    
    let frame = xyzt::ReferenceFrames::InertialCartesian;
    let dv_dpt = dv1 - orbit_i.calc_motion(start_date_time.second() as f64, frame)[1];
    let dv_arr = dv2 - orbit_f.calc_motion(stop_date_time.second() as f64, frame)[1];
    let dv_tot = dv_dpt.norm() + dv_arr.norm();
    return dv_tot
}


/// Izzo based lambert solver
///
/// Inputs
/// ------
/// grav_param: `decimal
///     Central Grav param of system
/// 
/// r_i: `[xyz,] decimal`
///     Initial position vector
/// 
/// r_f: `[xyz,] decimal`
///     Final position vector
/// 
/// tof: `decimal`
///     Time of flight in seconds
/// 
/// prograde_sign: `f64`
///     Prograde vs retrograde
/// 
/// Returns
/// -------
/// 
/// Sources
/// ------- 
/// link to polisatro's repo
/// 
/// link to izzo's paper
/// 
pub fn lambert_solve(
    grav_param: f64,
    r_i: Vector3<f64>,
    r_f: Vector3<f64>,
    tof: f64,
    prograde_sign: f64
) -> (Vector3<f64>, Vector3<f64>) {
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
    let mut lambda: f64 = prograde_sign * (1.0 - min).sqrt();
    let mut u_t_i: Vector3<f64>;
    let mut u_t_f: Vector3<f64>;

    // Orient tangential vectors
    if u_h[2] < 0.0 {
        lambda = -lambda;
        u_t_i = prograde_sign * u_r_i.cross(&u_h);
        u_t_f = prograde_sign * u_r_f.cross(&u_h);
    }
    else{
        u_t_i = prograde_sign * u_h.cross(&u_r_i);
        u_t_f = prograde_sign * u_h.cross(&u_r_f);
    }

    // Dimensionless time parameters
    let time = tof * ((2. * grav_param / semi_perim.powi(3)).sqrt());
    
    let x = find_xy(time, lambda);
    let y: f64 = (1. - lambda.powi(2) * (1. - x.powi(2))).sqrt();

    let gamma: f64 = (grav_param * semi_perim / 2.).sqrt();
    let rho: f64 = (r_i_norm - r_f_norm) / c_norm;
    let sigma: f64 = (1. - rho.powi(2)).sqrt();

    // Radial and Tangential Components
    let v_r_i = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / r_i;
    let v_r_f = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / r_f;
    let v_t_i = gamma * sigma * (y + lambda * x) / r_i;
    let v_t_f = gamma * sigma * (y + lambda * x) / r_f;

    // Velocity vectors
    let v_i: Vector3<f64> = v_r_i * u_r_i + v_t_i * u_t_i;
    let v_f: Vector3<f64> = v_r_f * u_r_f + v_t_f * u_t_f;

    return (v_i, v_f)
}

///  Find x and y in time lambda space (izzo 2014)
///
/// Inputs
/// ------
/// time: `decimal,  seconds`
/// 
/// lambda: `decimal, seconds`
/// 
fn find_xy(time: f64, lambda: f64) ->  Vector3<f64>{
    let mut M_max: f64 = time / PI; // floor to int
    let T_00: f64 = lambda.acos() + lambda * (1. - lambda.powi(2)).sqrt();

    // Refine maximum number of revolutions if necessary
    if (time < T_00 + M_max * PI) && (M_max > 0.0){
        
        //  Halley iteration
        //     _, T_min = _compute_T_min(ll, M_max, numiter, rtol)
        if time < T_min{M_max = M_max - 1.0;}
    }

    let T_1 = 2/3 * (1- lambda.powi(3));
    let mut x_0: f64;
    if time < T_1{
        x_0 = 2 * T_1/time  - 1; 
    }
    else if time >= T_1 &&  time < T_0{
        x_0 = (T_0  / time).powf((T_1/T_0).log2())
    } else {
        x_0 = T_0 / time - 1;
    }
    
    // FIXME-TD: Translate into Rust -V
    // Start Householder iterations from x_0 and find x, y
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
    //             psi = (x * y + ll * (1 - x**2)).acos();
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
}

#[cfg(test)]
mod interplan_tests {

    use super::*;
    use crate::donnager::constants as cst;

    #[test]
    fn test_calc_esc_vel()
    {
        let grav_param: f64 = cst::SUN::GRAV_PARAM;
        let orb_radius_0: f64 = 1000000.;
        let orb_radius_f: f64 = 2000000.;
        let esc_vel: f64 = calc_esc_vel(grav_param, orb_radius_0, orb_radius_f);
        assert_eq!(esc_vel, 13302555.410020536);
    }

    #[test]
    fn test_calc_mission_delta_v(){
        // let earth_orbit
        // let mars_orbit
        // let start_date
        // let stop_date
    }


}
