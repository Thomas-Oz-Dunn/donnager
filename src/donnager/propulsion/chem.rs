/*
Chemistry calculations
*/
use crate::donnager::propulsion::atom;

pub struct Chemical{
    name: String,
    moles: f64,
    formula: String,
    composition: Box<[(atom::Element, i32)]>,
    molar_mass: f64,
    heat_capacity_ratio: f64,
}

// pub fn calc_combustion_heat(
//     reactants: Box<[Chemical]>
// ) -> f64 {
//     // Vectorized calculation
//     // Minimize Gibbs free energy
// }

pub fn calc_chem_molar_mass(
    composition: Box<[(atom::Element, i32)]>
) -> f64 {
    let mut chem_molar_mass: f64 = 0.0;
    for i_element in composition.iter(){
        chem_molar_mass += i_element.0.molar_mass * f64::from(i_element.1);
    };
    return chem_molar_mass;
}