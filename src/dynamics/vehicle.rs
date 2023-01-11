/*
Vehicle modelling
*/

use crate::{propulsion as prop, cosmos};

#[derive(Clone, Debug, PartialEq)]
pub struct Vehicle{
    pub name: String,
    pub mass: f64,
    pub engine: prop::engine::Engine
}

#[derive(Clone, Debug, PartialEq)]
pub struct Multistage{
    pub name: String,
    pub stages: Vec<Vehicle>
}

impl Vehicle {
    /// Calculate fuel mass to reach a particular delta v
    pub fn calc_mass_fuel(
        &self, 
        delta_v: f64, 
        grav_acc: f64, 
        external_mass: f64
    ) -> f64 {
        let mass_ratio: f64 = prop::ballistics::calc_mass_ratio(
            delta_v, 
            self.engine.isp,
            grav_acc);
        let mass_fuel: f64 = (self.mass + external_mass) * mass_ratio;    
        return mass_fuel
    }
}


impl Multistage {

    /// Calculate total rest mass
    pub fn calc_inertial_mass(&self) -> f64 {
        let mut mass_mut: f64 = 0.0;
        self.stages.iter().for_each(|stage| mass_mut += stage.mass);
        let mass: f64 = mass_mut;
        return mass
    }


    /// Calculate total fuel mass required for delta v
    pub fn calc_mass_fuel(
        &self, 
        delta_v: f64, 
        launch_site: cosmos::space::SurfacePoint
    ) -> Vec<f64> {
        // Starting conditions
        let n_stage: usize = self.stages.len();
        let total_mass_0: f64 = self.calc_inertial_mass();
        let mut radius: f64 = launch_site.calc_surface_radius();
        let body: cosmos::grav::Body = launch_site.body;
        let mut grav_acc: f64 = body.calc_grav_acc(radius);

        let mut fuel_mass_mut: Vec<f64> = Vec::<f64>::new();
        for (stage, i_stage) in self.stages.iter().zip((0..n_stage).into_iter()) {
            
            let mut external_mass: f64 = 0.0;
            if i_stage > 0 {
                (0..i_stage).rev().for_each(|i_mass|
                    external_mass += self.stages[i_mass].mass + fuel_mass_mut[i_mass])
            };
            
            // For each stage, break the total delta v down proportionally
            let delta_v_stage: f64 = 
                delta_v * (stage.mass + external_mass) / (total_mass_0 + external_mass);
            fuel_mass_mut[i_stage] = stage.calc_mass_fuel(
                delta_v_stage, 
                grav_acc, 
                external_mass);
                
            // Calculate new starting conditions for next stage
            let mass_ratio: f64 = prop::ballistics::calc_mass_ratio(
                delta_v_stage, 
                stage.engine.isp, 
                grav_acc);
            radius = prop::ballistics::calc_burnout_height(
                mass_ratio, 
                grav_acc, 
                stage.engine.isp);
            grav_acc = body.calc_grav_acc(radius);
        }

        let fuel_mass: Vec<f64> = fuel_mass_mut;
        return fuel_mass;
    }

}