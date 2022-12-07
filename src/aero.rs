
#[path="./constants.rs"] mod constants;

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

pub fn calc_mach_number(
    velocity: f64,
    temperature: f64,
    gamma: f64
) -> f64 {
    let speed_of_sound: f64 = (gamma * constants::GAS_CONST * temperature).sqrt();
    let mach_number: f64 = velocity / speed_of_sound;
    mach_number 
}

