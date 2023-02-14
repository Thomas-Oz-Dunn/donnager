/*
Solar power generation simulation for sun-synchronous orbit
*/

use donnager::constants as cst;
use donnager::cosmos as cosm;
use donnager::dynamics as dynam;

let sample_tle: &str = 
    "ICEYE X1
    1 43114U 18004D   23045.46121586  .00040042 00000-0  11420-2 0  9993
    2 43114  97.3365 121.5406 0007153 204.3285 155.7622 15.35923585283273";
    
let sso_orbit: Orbit = cosm::gravity::Orbit.from_tle(sample_tle);

// Calculate occlusion events

let pos = sso_orbit.;
// Calculate flux over time

let solar_pressure = cosm::electromagnetism::calc_em_radiative_force()
// Replicate / Falsify the results Dae Lee's paper, cite?
