/*
Chemistry calculations
*/

pub struct Chemical{
    name: String,
    formula: String,
    composition: [tuple(String, ui32)],
    molar_mass: f64,
    thermal_conductivity: f64,
    heat_capacity_ratio: f64
}

pub fn calc_combustion_heat(
    reactants: [Chemical],
    products: [Chemical]
) -> f64 {
    // Vectorized calculation
}