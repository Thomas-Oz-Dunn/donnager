

use crate::propulsion as prop;

#[derive(Clone, Debug, PartialEq)]
pub struct Vehicle{
    pub name: String,
    pub mass_0: f64,
    pub mass_prop: f64,
    pub engine_type: String,
    pub engine_isp: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Multistage{
    pub name: String,
    pub stages: Vec<Vehicle>
}

impl Vehicle {
    /// Calculate fuel mass to reach delta v
    pub fn calc_mass_fuel(&self, delta_v: f64, grav_acc: f64) -> f64 {
        let engine_isp: f64 = self.engine_isp;
        let mass_ratio: f64 = prop::ballistics::calc_mass_ratio(delta_v, engine_isp, grav_acc);
        let mass_fuel: f64 = self.mass_0 * mass_ratio;    
        return mass_fuel
    }
}


impl Multistage {
    /// Calculate total fuel mass
    pub fn calc_mass_fuel(&self, delta_v: f64, grav_acc: f64) -> Vec<f64> {
        // for each stage, break the total delta v down

    }

}