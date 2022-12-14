/*
Chemistry calculations
*/

#[path="./atom.rs"] mod atom;

pub struct Chemical{
    name: String,
    moles: f64,
    formula: String,
    composition: [tuple(atom::Molecule, ui32)],
    molar_mass: f64,
    thermal_conductivity: f64,
    heat_capacity_ratio: f64,
    triple_point: tuple(f64, f64)
}

pub fn calc_combustion_heat(
    reactants: [Chemical]
) -> f64 {
    // Vectorized calculation
    // Minimize Gibbs free energy
}

pub fn mixture(
    mix: [Chemical]
) -> Chemical {
    // First, check reaction
    // If reaction, include products?
    let mut mix_name: String = "";
    for chemical in mix{
        mix_name += chemical.name + " ";
    }
    // Place all the mole counts and heat capacities in Vectors
    // Vectorized multiplication instead of serial

    let mixture = Chemical {
        name: mix_name,
    };
}