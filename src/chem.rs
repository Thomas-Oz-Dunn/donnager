/*
Chemistry calculations
*/

pub struct Chemical{
    name: String,
    moles: f64,
    formula: String,
    composition: [tuple(String, ui32)],
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
    Chemicals: [Chemical]
) -> Chemical {
    // First, check reaction
    // If reaction, include products?
    let mut mix_name: String = "";
    for chem in Chemicals{
        mix_name += chem.name + " ";
    }
    // Place all the mole counts and heat capacities in Vectors
    // Vectorized multiplication instead of serial

    let mixture = Chemical {
        name: mix_name,
    };
}