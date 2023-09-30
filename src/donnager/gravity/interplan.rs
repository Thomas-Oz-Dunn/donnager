/*
Interplanetary Mission Planner
*/
use std::f64::consts::PI;
use nalgebra::{Vector3};
use rayon::prelude::*;

use crate::donnager::{
    gravity::kepler as kepler, 
    spacetime as xyzt,
    math as math
};

/// Calculate interplanetary mission escape velocity
/// 
/// Inputs
/// ------
/// grav_param: `f64`
///     Gravitational parameter of centeral bodys
/// 
/// orb_radius_0: `f64`
///     Starting planet heliocentric radius
/// 
/// orb_radius_f: `f64`
///     Ending planet heliocentric radius
/// 
/// Outputs
/// -------
/// esc_vel: `f64`
///     Escape velocity
pub fn calc_esc_vel(
    grav_param: f64,
    orb_radius_0: f64,
    orb_radius_f: f64,
) -> f64 {
    let mean_radius: f64 = (orb_radius_0 + orb_radius_f) / 2.;
    return (grav_param * (2. / orb_radius_0 - 1. / mean_radius)).sqrt();
}

// /// Calculate delta v for mission in range
// /// 
// /// Inputs                                                   
// /// ------
// /// start_date_time: `DateTime<Utc>`
// ///     Start datetime for mission
// /// 
// /// stop_date_time: `DateTime<Utc>`
// ///     stop datetime for mission
// /// 
// /// orbit_1: `Orbit`
// ///     Orbit of starting planet
// /// 
// /// orbit_2: `Orbit`
// ///     Orbit of ending planet
// /// 
// /// n_evals: `u32`
// /// 
// /// rtol: `f64`
// /// 
// /// numiter: `i32`
// pub fn calc_mission_delta_v(
//     start_date_time: DateTime<Utc>,
//     stop_date_time: DateTime<Utc>,
//     orbit_i: kepler::Orbit,
//     orbit_f: kepler::Orbit,
//     n_evals: u32,
//     rtol: f64,
//     numiter: i32
// ) -> f64 {
//     let t_start: f64 = start_date_time.timestamp() as f64;
//     let t_end: f64 = stop_date_time.timestamp() as f64;
//     let time_step: f64 = (t_end - t_start) / n_evals as f64;
//     let eval_times: [f64; 2] = [t_start, t_end];  
//     for i_val in 0..n_evals{
//         let tof: f64 = time_step * i_val as f64;
//         let dvs: Vec<(Vector3<f64>, Vector3<f64>)> = parallel_lambert(
//             orbit_i, 
//             orbit_f, 
//             rtol, 
//             numiter, 
//             &eval_times,
//             tof
//         );
//         // Iterate across eval times
//         // 
//         // let dv_tot = dvs[i_time, 0].norm() + dvs[i_time, 1].norm();
//     }


//     return dv_tot
// }

pub fn parallel_lambert(
    orbit_i: kepler::Orbit, 
    orbit_f: kepler::Orbit, 
    rtol: f64, 
    numiter: i32, 
    eval_times:  &[f64],
    tof:  f64
) -> Vec<(Vector3<f64>, Vector3<f64>)> {
    // TODO-TD: parallel across orbits 
    let ref_frame: xyzt::ReferenceFrames = xyzt::ReferenceFrames::InertialCartesian;
    
    let pos_idx: i8 = 0;
    let vel_idx: i8 =  1;

    // Parallel across start time 
    let parallel_val_times = eval_times.par_iter()
        .map(|  eval_time  | {
        
        // TODO-TD: parallel across tof

        let r_i: Vector3<f64> = orbit_i.calc_motion(
            *eval_time, 
            ref_frame, 
            pos_idx
        );
        let v_i: Vector3<f64> = orbit_i.calc_motion(
            *eval_time, 
            ref_frame, 
            vel_idx
        );

        let r_f: Vector3<f64> = orbit_f.calc_motion(
            *eval_time + tof, 
            ref_frame, 
            pos_idx
        );
        let v_f: Vector3<f64> = orbit_f.calc_motion(
            *eval_time + tof, 
            ref_frame, 
            vel_idx
        );

        let helio_vs: (Vector3<f64>, Vector3<f64>) = lambert_solve(
            orbit_i.central_body.grav_param,  
            r_i,  
            r_f, 
            tof, 
            1.0,
            false,
            rtol,
            numiter
        );

        let dv_dpt: Vector3<f64> = helio_vs.0 - v_i;
        let dv_arr: Vector3<f64> = helio_vs.1 - v_f;
                    
        return (dv_dpt, dv_arr)
    });

    let dvs: Vec<_> = parallel_val_times.collect();
    return dvs
        
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
/// v_i: `[xyz,]  decimal`
///     Initial Velocity vector
/// 
/// v_f: `[xyz,]  decimal`
///     Final Velocity vector
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
    is_max: bool,
    rtol: f64,
    n_iter: i32
) -> (Vector3<f64>, Vector3<f64>) {
    // Constants
    let x_idx: usize = 0;
    let y_idx: usize = 1;
    let z_idx: usize = 2;

    // Normalize operations
    let chord: Vector3<f64> = r_f - r_i;
    let c_norm: f64 = chord.norm();
    let r_i_norm: f64 = r_i.norm();
    let r_f_norm: f64 = r_f.norm();
    let semi_perim: f64 = (r_i_norm + r_f_norm + c_norm) * 0.5;

    let u_r_i: Vector3<f64> = r_i / r_i_norm;
    let u_r_f: Vector3<f64> = r_f / r_f_norm;

    let i_h: Vector3<f64> = u_r_i.cross(&u_r_f);
    let u_h: Vector3<f64> = i_h / i_h.norm();
    
    // Orient tangential vectors
    let is_clockwise: bool = u_h[z_idx] < 0.0;
    let min: f64 = (1. as f64).min(c_norm / semi_perim);
    let mut lambda_mut: f64 = prograde_sign * (1.0 - min).sqrt();
    let u_t_i_mut: Vector3<f64>;
    let u_t_f_mut: Vector3<f64>;

    if is_clockwise {
        lambda_mut = -lambda_mut;
        u_t_i_mut = prograde_sign * u_r_i.cross(&u_h);
        u_t_f_mut = prograde_sign * u_r_f.cross(&u_h);

    } else {
        u_t_i_mut = prograde_sign * u_h.cross(&u_r_i);
        u_t_f_mut = prograde_sign * u_h.cross(&u_r_f);

    }

    let lambda: f64 = lambda_mut;
    let u_t_i: Vector3<f64> = u_t_i_mut;
    let u_t_f: Vector3<f64> = u_t_f_mut;

    // Dimensionless time parameters
    let time: f64 = tof * ((2. * grav_param / semi_perim.powi(3)).sqrt());
    let mean_motion: f64 = time / PI;
    let xy: Vec<f64> = find_xy(
        time, 
        mean_motion, 
        lambda, 
        is_max, 
        rtol, 
        n_iter
    ).unwrap();

    let x: f64 = xy[x_idx];
    let y: f64 = xy[y_idx];

    let gamma: f64 = (grav_param * semi_perim / 2.).sqrt();
    let rho: f64 = (r_i_norm - r_f_norm) / c_norm;
    let sigma: f64 = (1. - rho.powi(2)).sqrt();
    let l_1: f64 = lambda * y - x;
    let l_2: f64 = lambda * y + x;

    // Radial and Tangential Components
    let v_r_i: f64 = gamma * (l_1 - rho * l_2) / r_i_norm;
    let v_r_f: f64 = -gamma * (l_1 + rho * l_2) / r_f_norm;
    let v_t_i: f64 = gamma * sigma * (y + lambda * x) / r_i_norm;
    let v_t_f: f64 = gamma * sigma * (y + lambda * x) / r_f_norm;

    // Velocity vectors
    let v_i: Vector3<f64> = v_r_i * u_r_i + v_t_i * u_t_i;
    let v_f: Vector3<f64> = v_r_f * u_r_f + v_t_f * u_t_f;

    return (v_i, v_f)
}

/// Find x and y in time lambda space (Izzo 2014)
///
/// Inputs
/// ------
/// time: `f64,  seconds`
/// 
/// mean_motion: `f64`
/// 
/// lambda: `decimal, seconds`
/// 
/// is_max
/// 
/// rtol
/// 
/// n_iter
/// 
fn find_xy(
    time: f64, 
    mean_motion: f64, 
    lambda: f64, 
    is_max: bool,
    rtol: f64,
    n_iter: i32
) ->  Result<Vec<f64>, &'static str>{
    let mut m_max: f64 = time / PI; // floor to int
    let t_00: f64 = _lambda_to_t_0(lambda);

    let mut t_min: f64 = 0.0;
    // Refine maximum number of revolutions if necessary
    if (time < t_00 + m_max * PI) && (m_max > 0.0){

        _compute_t_min(
            lambda, 
            &mut t_min, 
            m_max, 
            rtol, 
            n_iter
        );

        if time < t_min{
            m_max = m_max - 1.0;
        }
    }

    if mean_motion > m_max{
        let _ = Err::<Vec<f64>, &str>("Mean motion > max M");
    }

    // Initial guess
    let x_0: f64 = _initial_guess(
        time, 
        lambda, 
        mean_motion, 
        is_max
    );

    // Start Householder iterations from x_0 and find x, y
    let x: f64 = householder(
        x_0, 
        time, 
        lambda, 
        mean_motion, 
        rtol, 
        n_iter
    );

    let y: f64 = _compute_y(x, lambda);
    Ok(vec![x, y])

}


/// Initial Guess for lambert solver
fn _initial_guess(
    time: f64, 
    lambda: f64, 
    mean_motion:f64, 
    is_max: bool
) -> f64 {
    let x_0: f64;

    if mean_motion == 0.0 {
        // Singular Revolution
        let t_00: f64 = _lambda_to_t_0(lambda);
        let t_0: f64 = t_00 + mean_motion * PI;  
        let t_1: f64 = (2. as f64).powf(1. - lambda.powi(3)) / 3.;

        if time >= t_0 {
            x_0 = (t_0 / time).powf(2. / 3.) - 1.;

        } else if time < t_1 {
            x_0 = 5. / 2. * t_1 / time * (t_1 - time) / (1. - lambda.powi(5)) + 1.;

        } else {
            // See https://github.com/poliastro/poliastro/issues/1362
            let base: f64 = (2. as f64).ln() * (time / t_0).ln() / (t_1 / t_0).ln();
            x_0 = (base).exp() - 1.;
        }
        
    } else {
        // Multiple revolutions
        let v_l: f64 = (
            ((mean_motion +  1.) * PI) / (8. * time)
        ).powf(2. / 3.);

        let v_r: f64 = (
            (8. * time) / (mean_motion * PI)
        ).powf(2. / 3.);   
        
        let x_0l: f64 = (v_l - 1.) / (v_l + 1.);
        let x_0r: f64 = (v_r - 1.) / (v_r + 1.);
        
        if is_max {
            x_0 = x_0l.max(x_0r);

        } else {
            x_0 = x_0l.min(x_0r);

        }
    }
    let x: f64 = x_0;
    return x        
}

/// Map between lambda and t_0
fn _lambda_to_t_0(lambda: f64) -> f64 {
    lambda.acos() + lambda * (1. - lambda.powi(2)).sqrt()
}

/// Computer minimum T
fn _compute_t_min(
    lambda: f64, 
    t_min: &mut f64, 
    m_max: f64, 
    rtol: f64,
    n_iter: i32
) {
    if lambda == 1.0{
        let x_t_min: f64 = 0.0;
        let y = _compute_y(x_t_min, lambda);
        *t_min = _tof_equation(
            x_t_min, 
            y,
            lambda, 
            m_max
        );
    } else {
        if m_max == 0.0 {
            *t_min = 0.0;
        }
        else{
            let x_i: f64 = 0.1;
            let y_i = _compute_y(x_i, lambda);

            let x_t_min: f64 = _halley(
                x_i, 
                y_i, 
                lambda, 
                m_max,
                rtol, 
                n_iter
            );

            let y_t_min = _compute_y(x_t_min, lambda);
            *t_min = _tof_equation(
                x_t_min, 
                y_t_min,
                lambda,
                m_max
            );
        }
    }
}

///
fn _halley(
    x_i: f64, 
    y_i: f64, 
    lambda: f64, 
    m_max: f64,
    rtol: f64, 
    n_iter: i32
) -> f64 {
    let mut x_iter = x_i;
    for _ in 1..n_iter
    {
        let t = _tof_equation(
            x_iter, 
            y_i,
            lambda,
            m_max
        );
        let fder = _tof_equation_d(
            x_iter, 
            y_i, 
            t, 
            lambda
        );
        let fder2 = _tof_equation_d2(
            x_iter, 
            y_i, 
            t, 
            fder, 
            lambda
        );

        let x: f64 = x_i - 2_f64 * t * fder / (2_f64 * (fder2).powi(2) - t * fder2);
        
        if (x - x_iter).abs() < rtol{
            return x
        }

        x_iter = x;
    }

    return x_iter

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
    y: f64,
    lambda: f64, 
    mean_motion: f64
)-> f64 {
    let time_: f64;

    let cond_0: bool = (0.6 as f64).sqrt() < x;
    let cond_1: bool = mean_motion == 0.0;
    let cond_2: bool = x < (1.4 as f64).sqrt();

    if cond_0 && cond_1 && cond_2{
        let eta: f64 = y - lambda * x;
        let s: f64 = (1. - lambda - x * eta) * 0.5;
        let q: f64 = 4. / 3. * math::hyp2f1b(s); 
        time_ = (eta.powi(3) * q + 4. * lambda * eta) * 0.5;

    } else {
        let psi: f64;
        if -1. <= x && x < 1. {
            psi = (x * y + lambda * (1. - x.powi(2))).acos();

        } else if x >= 1. {
            psi = ((y - x * lambda) * (x.powi(2) - 1.).sqrt()).asinh();

        } else {
            psi = 0.0;

        }
        
        let z: f64 = 1. - x.powi(2);
        let sqrt_z:  f64 = z.abs().sqrt();
        time_ = (((psi + mean_motion * PI)/sqrt_z) - x + lambda * y) / z
        
    }
    return time_

}

/// Householder iteration
fn householder(
    x_0: f64, 
    time: f64, 
    lambda: f64, 
    mean_motion: f64, 
    rtol: f64,
    n_iter: i32
) -> f64 {
    let mut x_0: f64 = x_0;

    for _ in 1..n_iter
    {
        let y: f64 = _compute_y(x_0, lambda);
        let fval = _tof_equation(
            x_0, 
            y,
            lambda, 
            mean_motion
        );

        let t: f64 = fval - time;
        let fder: f64 = _tof_equation_d(
            x_0, 
            y, 
            t, 
            lambda
        );
        let fder2: f64 = _tof_equation_d2(
            x_0, 
            y, 
            t, 
            fder, 
            lambda
        );
        let fder3: f64 = _tof_equation_d3(
            x_0, 
            y, 
            fder, 
            fder2, 
            lambda
        );

        // Householder step (quartic)
        let numer: f64 = fder.powi(2) - fval * fder2 / 2.;

        let d_1: f64 = fder * (fder.powi(2) - fval * fder2);
        let denom: f64 = d_1 + fder3 * fval.powi(2) / 6.;

        let x: f64 = x_0 - fval * numer / denom;

        if (x - x_0).abs() < rtol{
            return x
        }
        x_0 = x

    }
    return x_0

}

fn _compute_y(x: f64, lambda: f64) ->  f64 {
    return (1. - lambda.powi(2) * (1. - x.powi(2))).sqrt()
}


/// Time of flight first order derivative
fn _tof_equation_d(
    x:  f64, 
    y:  f64, 
    time:  f64, 
    lambda:  f64
) -> f64 {
    let numer: f64 = 3. * time * x - 2. + 2. * lambda.powi(3) * x / y;
    return  numer / (1. - x.powi(2));
}


/// Time of flight second order derivative
fn _tof_equation_d2(
    x_0:  f64, 
    y:  f64, 
    time:  f64, 
    fder: f64, 
    lambda:  f64
) -> f64{
    let a: f64 = 3. * time;
    let b: f64 = 5. * x_0 * fder;
    let c_1: f64 = (2. as f64).powf(1. - lambda.powi(2));
    let c: f64 = c_1 * lambda.powi(3) / y.powi(3);
    return (a + b + c) / (1. - x_0.powi(2))
}


/// Time of flight third order derivative
fn _tof_equation_d3(
    x_0: f64, 
    y: f64, 
    fder: f64, 
    fder2: f64, 
    lambda: f64
) -> f64{
    let a: f64 = 7. * x_0 * fder2;
    let b: f64 = 8. * fder;
    let c: f64 = 6. * (1. - lambda.powi(2)) * lambda.powi(5) * x_0 / y.powi(5);
    return  (a + b - c) / (1. - x_0.powi(2));
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
        let esc_vel: f64 = calc_esc_vel(
            grav_param, 
            orb_radius_0, 
            orb_radius_f
        );
        assert_eq!(esc_vel, 13302555.410020536);
    }

    // #[test]
    // fn test_calc_mission_delta_v(){
    //     // let earth_orbit
    //     // let mars_orbit
    //     // let start_date
    //     // let stop_date
    // }

    #[test]
    fn test_find_xy(){
        let time = 0.;
        let mean_motion = 50.234;
        let lambda = 0.;
        let is_max = true;
        let rtol = 0.00000001;
        let n_iter  = 1000;

        let xy = find_xy(
            time,
            mean_motion,
            lambda,
            is_max,
            rtol,
            n_iter
        );
    }


    #[test]
    fn test_lambert_solve(){

        let grav_param = cst::EARTH::GRAV_PARAM;
        let r_i: Vector3<f64> = Vector3::from_vec(vec![0. as f64, 0. as f64, 1.5 * cst::EARTH::RADIUS_EQUATOR]);
        let r_f: Vector3<f64> = Vector3::from_vec(vec![0. as f64, 1.4 * cst::EARTH::RADIUS_EQUATOR, 0. as f64]);
        let tof  = 235.8776;
        let prograde_sign = 1.;
        let is_max= true;
        let rtol = 1e-6;
        let n_iter = 1000;

        let solve = lambert_solve(
            grav_param,
            r_i,
            r_f,
            tof,
            prograde_sign,
            is_max,
            rtol,
            n_iter,
        );
    }
}
