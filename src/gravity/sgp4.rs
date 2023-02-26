/* 
Propogator scrap file
*/
use crate::gravity::kepler as kepler;

// State
// Objects- pos, vel, mass, grav field model

pub struct GravConstants {
    pub geopotential: Geopotential,
    pub raan_dot: f64,
    pub argument_of_perigee_dot: f64,
    pub mean_anomaly_dot: f64,
    pub c1: f64,
    pub c4: f64,
    pub k0: f64,
    pub k1: f64,
    pub orbit_0: kepler::Orbit,
}


pub struct Geopotential {
    pub eq_radius: f64, // Equatorial radius
    pub ke: f64, // square root of earth's grav param in earth radii^3 min^-2
    pub j2: f64, // un-normalised second zonal harmonic
    pub j3: f64, // un-normalised third zonal harmonic
    pub j4: f64, // un-normalised fourth zonal harmonic
}

pub const WGS84: Geopotential = Geopotential {
    eq_radius: 6378.137,
    ke: 0.07436685316871385,
    j2: 0.00108262998905,
    j3: -0.00000253215306,
    j4: -0.00000161098761,
};

