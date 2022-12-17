/*
Ballistics calculations
*/

// TODO-TD: mutlistage trade studies

pub fn calc_mass_ratio(
    delta_v: f64,
    engine_isp: f64,
    grav_acc: f64
) -> f64{
    let v_exhaust: f64 = engine_isp * grav_acc;
    let power: f64 = delta_v / v_exhaust;
    let mass_ratio: f64 = power.exp() - 1.0;
    mass_ratio
}

pub fn calc_burnout_height(
    mass_ratio: f64,
    grav_acc: f64,
    engine_isp: f64
) -> f64{
    let z: f64 = 1.0 - mass_ratio;
    let y: f64 = z * (z).ln() + mass_ratio - mass_ratio.powi(2) / 2.0;
    let burnout_height: f64 = grav_acc * engine_isp.powi(2) * y;
    burnout_height
}


pub fn calc_coast_height(
    mass_ratio: f64,
    acc_ratio: f64,
    engine_isp: f64,
    grav_acc: f64
) -> f64 {
    let z: f64 = 1.0 - mass_ratio;
    let x: f64 = (1.0 / (z)).ln() - 1.0 / acc_ratio;
    let v_bo: f64 = grav_acc * engine_isp * (x);
    let coast_height: f64 = v_bo.powi(2) / (2.0 * grav_acc);
    coast_height
}
