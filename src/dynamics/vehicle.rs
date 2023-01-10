

use crate::{propulsion::{ballistics, engine}, cosmos::space::SurfacePoint};

#[derive(Clone, Debug, PartialEq)]
pub struct Vehicle{
    pub name: String,
    pub mass: f64,
    pub engine: engine::Engine
}

#[derive(Clone, Debug, PartialEq)]
pub struct Multistage{
    pub name: String,
    pub stages: Vec<Vehicle>
}

impl Vehicle {
    /// Calculate fuel mass to reach delta v
    pub fn calc_mass_fuel(&self, delta_v: f64, grav_acc: f64, external_mass: f64) -> f64 {
        let engine_isp: f64 = self.engine.isp;
        let mass_ratio: f64 = ballistics::calc_mass_ratio(delta_v, engine_isp, grav_acc);
        let mass_fuel: f64 = (self.mass + external_mass) * mass_ratio;    
        return mass_fuel
    }
}


impl Multistage {
    /// Calculate total fuel mass
    pub fn calc_mass_fuel(&self, delta_v: f64, launch_site: SurfacePoint) -> Vec<f64> {
        let n_stage: usize = self.stages.len();
        let mut radius: f64 = launch_site.calc_surface_radius();
        let mut grav_acc: f64 = launch_site.body.calc_grav_acc(radius);
        let mut fuel_mass_mut: Vec<f64> = Vec::<f64>::new();
        let mut total_mass_0: f64 = 0.0;

        // Calculate total mass
        for stage in self.stages.iter().rev() {
            total_mass_0 += stage.mass;
        }

        for (stage, i_stage) in self.stages.iter().zip((0..n_stage).into_iter()) {
            
            let mut external_mass: f64 = 0.0;
            if i_stage > 0 {
                for i_mass in (0..i_stage).rev(){
                    external_mass += self.stages[i_mass].mass + fuel_mass_mut[i_mass];
                }
            };
            
            // For each stage, break the total delta v down proportionally
            let delta_v_stage: f64 = delta_v * (stage.mass + external_mass) / total_mass_0;
            let mass_ratio: f64 = ballistics::calc_mass_ratio(
                delta_v_stage, 
                stage.engine.isp, 
                grav_acc);

            radius = ballistics::calc_burnout_height(mass_ratio, grav_acc, stage.engine.isp);
            grav_acc = launch_site.body.calc_grav_acc(radius);
            fuel_mass_mut[i_stage] = stage.calc_mass_fuel(
                delta_v_stage, 
                grav_acc, 
                external_mass);
        }

        let fuel_mass: Vec<f64> = fuel_mass_mut;
        return fuel_mass;
    }

}