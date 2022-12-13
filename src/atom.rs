/*
Atomic physics functions
*/

#[path="./constants.rs"] mod constants;

pub struct Element{
    name: String,
    number: ui32,
    molar_mass: f64,
    triple_point: tuple(f64, f64),
    critical_point: tuple(f64, f64),
    half_life: f64,
    van_der_waals_radius: f64,
    ionization_energies: [f64],
}

pub struct Molecule{
    element: Element,
    isotope: ui32,
    phase: String,
    decay_type: String
}

pub fn calc_binding_energy(
    mol: Molecule
) -> f64 {
    let num_protons: ui32 = mol.element.number;
    let num_neutrons: ui32 = mol.isotope = num_protons;
    let mass_protons: f64 = f64::from(num_protons) * constants::MASS_PROTON;
    let mass_protons: f64 = f64::from(num_neutrons) * constants::MASS_NEUTRON;
    let mass_defect: f64 = mass_molecule - (mass_protons + mass_neutrons);
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