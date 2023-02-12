/*
Solar power generation simulation for sun-synchronous orbit
*/

use donnager::constants as cst;
use donnager::cosmos as cosm;
use donnager::dynamics as dynam;


// TODO-TD: find TLE for sample SSO satellites
let sample_tle = '';
let sso_orbit: Orbit = cosm::gravity::Orbit.from_tle(sample_tle);

// Calculate occlusion events

// Calculate flux over time

// Replicate / Falsify the results Dae Lee's paper, cite?
