/*
Chemistry calculations
*/

#[path="./atom.rs"] mod atom;

pub struct Chemical{
    name: String,
    moles: f64,
    formula: String,
    composition: [(atom::Element, i32)],
    molar_mass: f64,
    heat_capacity_ratio: f64,
}

pub fn calc_combustion_heat(
    reactants: [Chemical]
) -> f64 {
    // Vectorized calculation
    // Minimize Gibbs free energy
}

pub fn calc_chem_molar_mass(
    composition: [(atom::Element, i32)]
) -> f64 {
    let mut chem_molar_mass: f64 = 0;
    for i_element in composition{
        chem_molar_mass += i_element[0].molar_mass * f64::from(i_element[1]);
    };
    return chem_molar_mass;
}