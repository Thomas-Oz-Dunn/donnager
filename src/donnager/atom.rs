/*
Atomic physics functions
*/

#[path="./constants.rs"] mod constants;

pub struct Element{
    pub name: String,
    pub number: i32,
    pub molar_mass: f64,
    pub isotope: i32,
    pub decay_type: String
}

pub fn calc_binding_energy(
    element: Element
) -> f64 {
    let num_protons: i32 = element.number;
    let num_neutrons: i32 = element.isotope - num_protons;
    let mass_protons: f64 = f64::from(num_protons) * constants::MASS_PROTON;
    let mass_neutrons: f64 = f64::from(num_neutrons) * constants::MASS_NEUTRON;
    let mass_defect: f64 = element.molar_mass - (mass_protons + mass_neutrons);
    let binding_energy: f64 = mass_defect * constants::SPEED_OF_LIGHT.powi(2);
    return binding_energy
}

pub fn rest_mass_to_ev(
    mass: f64
) -> f64 {
    let rest_energy: f64 = mass * constants::SPEED_OF_LIGHT.powi(2);
    let rest_ev: f64 = rest_energy * constants::ELECTRON_VOLT;
    return rest_ev
}

// TODO-TD: autogen every element and isotope from a database
pub const HYDROGEN_1: Element = Element {
    name: *"Hydrogen",
    number: 1,
    molar_mass: 1.0,
    isotope: 1,
    decay_type: *"stable"
};

pub const DEUTRIUM: Element = Element {
    name: *"Deutrium",
    number: 1,
    molar_mass: 2.014,
    isotope: 2,
    decay_type: *"stable"
};

pub const OXYGEN: Element = Element {
    name: *"Oxygen",
    number: 8,
    molar_mass: 16.0,
    isotope: 16,
    decay_type: *"stable"
};