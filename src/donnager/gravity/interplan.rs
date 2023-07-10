/*
Interplanetary Mission Planner
*/
use std::{f64::consts::PI, error::Error};
use chrono::{DateTime, Utc, Timelike};
use nalgebra::{Vector3};

use crate::donnager::{
    gravity::kepler as kepler, 
    spacetime as xyzt,
    math as math};

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
    let prograde_sign = 1.0;
    let is_max = false;

    let (dv1, dv2) = lambert_solve(
        grav_param,  
        r_i,  
        r_f, 
        tof, 
        prograde_sign,
        is_max
    );
        
    
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
/// v_i: `[xyz, ]  decimal`
/// 
/// v_f: `[xyz, ]  decimal`
/// 
/// Sources
/// ------- 
/// Poliastro (Python Implementation)
/// https://github.com/poliastro/poliastro
/// 
/// Dario Izzo, Revisiting Lambert's Problem. 2014 
/// https://arxiv.org/abs/1403.2705
/// 
pub fn lambert_solve(
    grav_param: f64,
    r_i: Vector3<f64>,
    r_f: Vector3<f64>,
    tof: f64,
    prograde_sign: f64,
    is_max: bool
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
    let time: f64 = tof * ((2. * grav_param / semi_perim.powi(3)).sqrt());
    let mean_motion = time / PI;
    let xy = find_xy(time, mean_motion, lambda,  is_max);
    let x = xy[0];
    let y = xy[1];
    let y: f64 = (1. - lambda.powi(2) * (1. - x.powi(2))).sqrt();

    let gamma: f64 = (grav_param * semi_perim / 2.).sqrt();
    let rho: f64 = (r_i_norm - r_f_norm) / c_norm;
    let sigma: f64 = (1. - rho.powi(2)).sqrt();

    // Radial and Tangential Components
    let v_r_i: f64 = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / r_i_norm;
    let v_r_f: f64 = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / r_f_norm;
    let v_t_i: f64 = gamma * sigma * (y + lambda * x) / r_i_norm;
    let v_t_f: f64 = gamma * sigma * (y + lambda * x) / r_f_norm;

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
fn find_xy(
    time: f64, 
    mean_motion: f64, 
    lambda: f64, 
    is_max: bool
) ->  Vec<f64>{
    let mut M_max: f64 = time / PI; // floor to int
    let T_00: f64 = lambda.acos() + lambda * (1. - lambda.powi(2)).sqrt();

    let mut T_min:  f64;
    // Refine maximum number of revolutions if necessary
    if (time < T_00 + M_max * PI) && (M_max > 0.0){
        _compute_t_min(lambda, &mut T_min, M_max);
        if time < T_min{
            M_max = M_max - 1.0;
        }
    }

    if mean_motion > M_max{
        // error out
    }

    //  Initial guess
    let x_0: f64 = _initial_guess(time, lambda, mean_motion, is_max);

    //  Start Householder iterations from x_0 and find x, y
    let x: f64 = householder(x_0, time, lambda, mean_motion, rtol, numiter);
    let y: f64 = _compute_y(x, lambda);
    return vec![x, y]

}

fn _initial_guess(
    time: f64, 
    lambda: f64, 
    mean_motion:f64, 
    is_max: bool
) -> f64 {
    let mut x_0: f64;
    if mean_motion == 0.0{
        // Singular Revolution
        let T_0: f64 = lambda.acos() + lambda * (1. - lambda.powi(2)).sqrt() + mean_motion * PI;  
        let T_1: f64 = (2. as f64).powf(1. - lambda.powi(3)) / 3.;
        if time >= T_0{
            x_0 = (T_0 / time).powf(2. / 3.) - 1.;
        } else if time < T_1{
            x_0 = 5. / 2. * T_1 / time * (T_1 - time) / (1. - lambda.powi(5)) + 1.;
        } else{
            // See https://github.com/poliastro/poliastro/issues/1362
            x_0 = ((2. as f64).ln() * (time / T_0).ln() / (T_1 / T_0).ln()).exp() - 1.;
        }
        
    } else  {
        // Multiple revolutions
        let v_l: f64 = ((mean_motion * PI + PI) / (8. * time)).powf(2. / 3.);
        let x_0l: f64 = (v_l - 1.) / (v_l + 1.);
        let v_r: f64= ((8. * time) / (mean_motion * PI)).powf(2. / 3.);   
        let x_0r: f64 = (v_r - 1.) / (v_r + 1.);
        
        if is_max{
            x_0 = x_0l.max(x_0r);
        }  else{
            x_0 = x_0l.min(x_0r);
        }
    }
    return x_0        
}

fn _compute_t_min(lambda: f64, T_min: &mut f64, M_max: f64) {
    if lambda == 1.0{
        let x_T_min = 0.0;
        *T_min = _tof_equation(
            x_T_min, 
            0.0, 
            lambda, 
            M_max);
    }else{
        if M_max == 0.0 {
            *T_min = 0.0;
        }
        else{
            let x_i: f64 = 0.1;
            let T_i: f64 = _tof_equation(
                x_i, 
                0.0, 
                lambda,
                M_max);
            let x_T_min = _halley(
                x_i, 
                T_i, 
                lambda, 
                rtol, 
                numiter);
            *T_min = _tof_equation(
                x_T_min, 
                0.0, 
                lambda,
                M_max);
        }
    }
}

/// Calculate time of flight equations
/// 
/// Parameters
/// ----------
/// x
/// 
/// T_0
/// 
/// lambda
/// 
/// mean_motion
/// 
/// Returns
/// -------
/// tof: `f64`
///     Time of flight
fn _tof_equation(
    x: f64, 
    T0: f64, 
    lambda: f64, 
    mean_motion: f64
)-> f64 {
    let y: f64 = (1. - lambda.powi(2) * (1. - x.powi(2))).sqrt();
    let mut time_: f64;

    if mean_motion == 0.0 && (0.6_f64).sqrt() < x && x < (1.4_f64).sqrt(){
        let eta:   f64 = y - lambda * x;
        let s: f64 = (1. - lambda - x * eta) * 0.5;
        let q:   f64 = 4. / 3. * math::hyp2f1b(s); 
        time_ = (eta.powi(3) * q + 4. * lambda * eta) * 0.5;

    } else {
        let mut psi: f64;
        if -1. <= x && x < 1.{
            psi = (x * y + lambda * (1. - x.powi(2))).acos();
        }
        else if x >= 1.{
            psi = ((y - x * lambda) * (x.powi(2) - 1.).sqrt()).asinh();
        }
        else{
            psi = 0.0;
        }
        let z: f64 = 1. - x.powi(2);
        let sqrt_z:  f64 = ((z).abs()).sqrt();
        time_ = (((psi + mean_motion * PI)/sqrt_z) - x + lambda * y) / z
        
    }
    return time_ - T0

}

/// Householder iteration
fn householder(
    x_0: f64, 
    time: f64, 
    lambda: f64, 
    mean_motion: f64, 
    rtol: f64,
    numiter: i32
) -> f64 {
    let mut x_0: f64  = x_0;
    for _ in 1..numiter
    {
        let y = _compute_y(x_0, lambda);
        let fval = _tof_equation_y(x_0, y, T0, ll, M);
        let T = fval + T0;
        let fder: f64 = _tof_equation_p(x_0, y, time, lambda);
        let fder2: f64 = _tof_equation_p2(x_0, y, time, fder, lambda);
        let fder3: f64 = _tof_equation_p3(x_0, y, time, fder, fder2, lambda);

        // Householder step (quartic)
        let x = x_0 - fval * (
            (fder.powi(2) - fval * fder2 / 2.)
            / (fder * (fder.powi(2) - fval * fder2) + fder3 * fval.powi(2) / 6.)
        );

        if (x - x_0).abs() < rtol{
            return x
        }
        x_0 = x
    }
    return Error

}

fn _compute_y(x: f64, lambda: f64) ->  f64 {
    return (1. - lambda.powi(2) * (1. - x.powi(2))).powf(0.5)
}

fn _tof_equation_p(x:  f64, y:  f64, time:  f64, lambda:  f64) -> f64 {
    return (3. * time * x - 2. + 2. * lambda.powi(3) * x / y) / (1. - x.powi(2));
}

fn _tof_equation_p2(x_0:  f64, y:  f64, time:  f64, fder: f64, lambda:  f64) -> f64{
    return (3. * time + 5. * x_0 * fder + (2. as f64).powf(1. - lambda.powi(2)) * lambda.powi(3) / y.powi(3)) / (1. - x_0.powi(2))
}

fn _tof_equation_p3(x_0: f64, y: f64, time: f64, fder: f64, fder2: f64, lambda: f64) -> f64{
    return (7. * x_0 * fder2 + 8. * fder - 6. * (1. - lambda.powi(2)) * lambda.powi(5) * x_0 / y.powi(5)) / (1. - x_0.powi(2));
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

    // #[test]
    // fn test_calc_mission_delta_v(){
    //     // let earth_orbit
    //     // let mars_orbit
    //     // let start_date
    //     // let stop_date
    // }


}
