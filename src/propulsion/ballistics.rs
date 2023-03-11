/*
Ballistics calculations
*/

// TODO-TD: mutlistage trade studies

/// Calculate mass ratio of payload to fuel
/// 
/// Inputs
/// ------
/// delta_v: `f64`
///     Delta-V of the payload in m/s.
/// 
pub fn calc_mass_ratio(
    delta_v: f64,
    isp: f64,
    grav_acc: f64
) -> f64{
    let v_exhaust: f64 = isp * grav_acc;
    let power: f64 = delta_v / v_exhaust;
    let mass_ratio: f64 = power.exp() - 1.0;
    mass_ratio
}

/// Calculate the burnout height of a rocket
pub fn calc_burnout_height(
    mass_ratio: f64,
    grav_acc: f64,
    isp: f64
) -> f64 {
    let z: f64 = 1.0 - mass_ratio;
    let y: f64 = z * (z).ln() + mass_ratio - mass_ratio.powi(2) / 2.0;
    let burnout_height: f64 = grav_acc * isp.powi(2) * y;
    burnout_height
}

/// Calculate coasting height of a rocket after engine cut off
pub fn calc_coast_height(
    mass_ratio: f64,
    acc_ratio: f64,
    isp: f64,
    grav_acc: f64
) -> f64 {
    let z: f64 = 1.0 - mass_ratio;
    let x: f64 = (1.0 / (z)).ln() - 1.0 / acc_ratio;
    let v_bo: f64 = grav_acc * isp * (x);
    let coast_height: f64 = v_bo.powi(2) / (2.0 * grav_acc);
    coast_height
}


#[cfg(test)]
mod ballistics_tests {
    use super::*;

    #[test]
    fn test_burnout(){
        let isp: f64 = 300.;
        let grav_acc: f64 = 9.81;
        let delta_v: f64 = 100.;
        let mass_ratio: f64 = calc_mass_ratio(
            delta_v, 
            isp, 
            grav_acc);
        let burnout_h: f64 = calc_burnout_height(
            mass_ratio, 
            grav_acc, 
            isp);
        assert_eq!(burnout_h, 6.1827991953331765);
    }
}
