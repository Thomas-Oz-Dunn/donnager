/*
Astrodynamics in Rust

*/
mod constants;

fn main() {
    let mass_0: f64 = constants::EARTH_MASS;
    let radius_0: f64 = constants::EARTH_RADIUS_EQUATOR; 

    let mass_1: f64 = 1.0; // kg
    let altitude: f64 = 3.0e6; // 300km, LEO
    let radius_f: f64 = radius_0 + altitude; // m
    let engine_isp: f64 = 300.0; // s

    let delta_v: f64 = orbital_velocity(mass_0, radius_f);
    let grav_acc: f64 = grav_acc(mass_0, radius_0);
    let mass_fuel: f64 = mass_fuel(delta_v, mass_1, engine_isp, grav_acc);    
    println!("{} kg of fuel", mass_fuel);
}

/*
TODO-TD: 
N-stage trade study
TLE ingest
Hohman Transfers
Trans Lunar Injections
Multithreading, Cloud Compute?
RK45 propogator
J2 perturbation
Launch cost calculator
Solar System Mineralogical data base query
Interplanetary Mission Plan (optimal launch windows, porkchop plot)
$ / kg (mineral X)
*/


pub fn grav_acc(
    mass_0: f64,
    radius: f64
) -> f64{
    let grav_acc: f64 = constants::GRAV_CONST * mass_0 / radius.powi(2);
    grav_acc
}

pub fn orbital_velocity(
    mass_0: f64,
    radius: f64
) -> f64{
    let energy: f64 = 2.0 * constants::GRAV_CONST * mass_0 / radius;
    let velocity: f64 = energy.sqrt();
    velocity
}

pub fn mass_fuel(
    delta_v: f64,
    mass_1: f64,
    engine_isp: f64,
    grav_acc: f64
) -> f64{
    let v_exhaust: f64 = engine_isp * grav_acc;
    let power: f64 = delta_v / v_exhaust;
    let mass_fuel: f64 = mass_1 * (power.exp() - 1.0);
    mass_fuel
}

pub fn apogee_height(
    mass_ratio: f64,
    engine_isp: f64,
    grav_acc: f64,
    acc_ratio: f64
) -> f64{
    let z: f64 = 1.0 - mass_ratio;
    let x: f64 = (1.0 / (z)).ln() - 1.0 / acc_ratio;
    let v_bo: f64 = grav_acc * engine_isp * (x);
    let coast_height: f64 = v_bo.powi(2) / (2.0 * grav_acc);
    
    let y: f64 = z * (z).ln() + mass_ratio - mass_ratio.powi(2) / 2.0;
    let burnout_height: f64 = grav_acc * engine_isp.powi(2) * y;
    let apogee_height: f64 = burnout_height + coast_height;
    apogee_height
}


pub fn calc_mass_flow(
    throat_area: f64,
    pressure_chamber: f64,
    temp_chamber: f64,
    gamma: f64,
    mol_weight: f64
) -> f64 {
    let a: f64 = (2.0/(gamma + 1.0)).powf((gamma + 1.0)/(2.0*(gamma - 1.0))); 
    let b: f64 = ((constants::GAS_CONST * temp_chamber) / gamma * mol_weight).sqrt();
    let c_star: f64 = b * a;

    let mass_flow: f64 = pressure_chamber * throat_area / c_star; 
    mass_flow
}