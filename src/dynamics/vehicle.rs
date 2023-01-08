

use crate::{propulsion::ballistics, cosmos::space::SurfacePoint};

#[derive(Clone, Debug, PartialEq)]
pub struct Vehicle{
    pub name: String,
    pub mass_0: f64,
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
        let mass_ratio: f64 = ballistics::calc_mass_ratio(delta_v, engine_isp, grav_acc);
        let mass_fuel: f64 = self.mass_0 * mass_ratio;    
        return mass_fuel
    }
}


impl Multistage {
    /// Calculate total fuel mass
    pub fn calc_mass_fuel(&self, delta_v: f64, launch_site: SurfacePoint) -> Vec<f64> {
        let mut radius_0: f64 = launch_site.calc_surface_radius();
        let mut grav_acc: f64 = launch_site.body.calc_grav_acc(radius_0);
        let mut fuel_mass_mut: Vec<f64> = Vec::<f64>::new();
        let mut i_stage: usize = 0;
        let mut total_mass_0: f64 = 0.0;

        for stage in self.stages.iter() {
            total_mass_0 += stage.mass_0;
        }

        for stage in self.stages.iter(){
            // for each stage, break the total delta v down proportionally
            let delta_v_stage: f64 = delta_v * stage.mass_0 / total_mass_0;
            let mass_ratio: f64 = ballistics::calc_mass_ratio(delta_v_stage, stage.engine_isp, grav_acc);
            let radius_f: f64 = ballistics::calc_burnout_height(mass_ratio, grav_acc, stage.engine_isp);
            
            grav_acc = launch_site.body.calc_grav_acc(radius_f);
            fuel_mass_mut[i_stage] = stage.calc_mass_fuel(delta_v_stage, grav_acc);
            radius_0 = radius_f;
            i_stage += 1;
        }

        let fuel_mass: Vec<f64> = fuel_mass_mut;
        return fuel_mass;
    }

}