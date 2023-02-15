/*
Solar power generation simulation for sun-synchronous orbit
*/
use nalgebra as na;
use na::Vector3;

use donnager::constants as cst;
use donnager::cosmos as cosm;

fn main(){
    let sample_tle: &str = 
        "ICEYE X1
        1 43114U 18004D   23045.46121586  .00040042 00000-0  11420-2 0  9993
        2 43114  97.3365 121.5406 0007153 204.3285 155.7622 15.35923585283273";
    
    let sso_orbit: cosm::gravity::Orbit = cosm::gravity::Orbit::from_tle(sample_tle.to_string());
    
    // Calculate occlusion events
    // ECI coordinates?
    // Where is sun in ECI @ ref datetime?
    
    let veh_pos_ecef = Vector3::<f64>.ones();
    let sun_pos_ecef = Vector3::<f64>.zeros();

    // Space params
    let src_obj_vec: Vector3<f64> = veh_pos_ecef - sun_pos_ecef;
    let em_flux: f64 = cst::SUN_MEAN_SOLAR_FLUX;

    // Vehicle paramsss
    let obj_reflecitivity: f64 = 1.6;
    let normal_area: f64 = 1;

    let solar_pressure = cosm::electromagnetism::calc_em_radiative_force(
        em_flux,
        src_obj_vec,
        obj_reflecitivity,
        normal_area
    );


    // Replicate / Falsify the results Dae Lee's paper, cite?

};
