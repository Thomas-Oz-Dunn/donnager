/*
Constants
*/

// TODO-TD: Provide sources for every constant

use std::f64::consts::PI;

/// Math
/// ----

/// Degrees to radians
pub const DEG_TO_RAD: f64 = PI / 180.0;
/// Radians to degrees
pub const RAD_TO_DEG: f64 = 180.0 / PI;

pub fn to_radians(degrees: f64) -> f64 {DEG_TO_RAD * degrees}
pub fn to_degrees(radians: f64) -> f64 {RAD_TO_DEG * radians}

pub const TERA: f64 = 1e12;
pub const GIGA: f64 = 1e9;
pub const MEGA: f64 = 1e6;
pub const KILO: f64 = 1e3;
pub const DECA: f64 = 1e1;

pub const DECI: f64 = 1e-1;
pub const CENTI: f64 = 1e-2;
pub const MILLI: f64 = 1e-3;
pub const MICRO: f64 = 1e-6;
pub const NANO: f64 = 1e-9;

/// Gravity
/// -------

/// Newtonian Gravitational constant
pub const GRAV_CONST: f64 = 6.6743015e-11; // m^3 * kg^-1 * s^-2

pub const J2000_DAY: f64 = 2451545.0;

pub struct SUN;
impl SUN {
    /// Mass of the Sun in kilograms
    pub const MASS: f64 = 1.9891e30; // kg
    /// Radius of the Sun in meters
    pub const RADIUS_EQUATOR: f64 = 9.6e7; // m
    pub const RADIUS_POLE: f64 = 9.5e7; // m
    /// Sun centered gravitational parameter
    pub const GRAV_PARAM: f64 = SUN::MASS * GRAV_CONST;
    pub const LUMINOSITY: f64 = 3.828e26; // Watts
    /// Mean Solar Flux at 1 AU
    pub const MEAN_SOLAR_FLUX: f64 = 1366.91; // Watts * m-2
    pub const ECC: f64 = 0.016709;
}

pub struct EARTH;
impl EARTH {
    /// Mass of the Earth in kilograms
    pub const MASS: f64 = 5.972e24; // kg
    /// Radius of Earth at equator in meters
    pub const RADIUS_EQUATOR: f64 = 6.378137e6; // m
    pub const RADIUS_POLE: f64 = 6.3567e6;
    
    /// Earth centered gravitational parameter
    pub const GRAV_PARAM: f64 = EARTH::MASS * GRAV_CONST; // m^3 * s^-2
    pub const ECC: f64 = 0.08182;
    pub const SIDEREAL_DAY: f64 = 23.9344696; // hours
    pub const ROT_RATE: f64 = 7.2921150e-5; // radians per second;
    pub const DAYS_PER_YEAR: f64 = 365.25;
    pub const AXIAL_TILT: f64 = -23.44;

}

/// Earth-Sun System
pub struct EarthSunOrbit;
impl EarthSunOrbit{
    pub const SEMI_MAJOR: f64 = 149.60e6; // km
    pub const ECC: f64 = 0.0167086;
    pub const ARG_PERIHELION: f64 = 102.9372;
    pub const INC: f64 = 0.0;
    pub const RAAN: f64 = 0.0;
    pub const MEAN_MOTION: f64 = 0.98560028; // radians per day
    pub const MEAN_ANOMALY: f64 = 0.0; // radians
}

pub struct MOON;
impl MOON{
    pub const MASS: f64 = 7.34767309e22; // kg
    pub const RADIUS_EQUATOR: f64 = 1.737e6; // m
    pub const GRAV_PARAM: f64 = MOON::MASS * GRAV_CONST; // m^3 * s^-2
}

// Moon-Earth System
pub struct MoonEarthOrbit;
impl MoonEarthOrbit{
    pub const ECC: f64 = 0.0549;
    pub const SEMI_MAJOR: f64 = 384400.0; // km   
    pub const INC: f64 = 5.145;
    pub const RAAN: f64 = 125.08;
    pub const ARG_PERIHELION: f64 = 318.15;
}


pub struct MARS;
impl MARS{
    pub const MASS: f64 = 6.4171e23; // kg
    pub const RADIUS_EQUATOR: f64 = 3.3895e6; // m
    pub const GRAV_PARAM: f64 = MARS::MASS * GRAV_CONST; // m^3 * s^-2
    pub const ECC: f64 = 0.0934;
    pub const SIDEREAL_DAY: f64 = 24.6229; // hours
    pub const ROT_RATE: f64 = 7.088e-5; // radians per second
}

// Mars-Sun System
pub struct MarsSunOrbit;
impl MarsSunOrbit{
    pub const SEMI_MAJOR: f64 = 227.9e6; // km
    pub const ECC: f64 = 0.0934;
    pub const ARG_PERIHELION: f64 = 286.5016;
    pub const INC: f64 = 1.850;
    pub const RAAN: f64 = 49.57854;
    pub const MEAN_MOTION: f64 = 0.524; // radians per day
    pub const MEAN_ANOMALY: f64 = 0.0; // radians
}


/// Electromagnetism
/// ----------------

/// Speed of light in a vacuum `m * s^-1`
pub const SPEED_OF_LIGHT: f64 = 2.99792458e8;
pub const ELECTRIC_CONST: f64 = 8.8541878123; // F * m^-1
pub const PLANCKS_CONST: f64 = 6.62607015e-34; // J * Hz-1
pub const BOLTZMAN_CONST: f64 = 1.380649e-23; // J * K-1

/// Nuclear Physics
/// ---------------

pub const MASS_PROTON: f64 = 1.67262e-27; //kg
pub const MASS_NEUTRON: f64 = 1.67492749e-27; //kg
pub const ELECTRON_VOLT: f64 = 1.602e-19; // Joule

/// Chemistry
/// ---------

pub const MOLE: f64 = 6.02e23; // molecules / mole
pub const GAS_CONST: f64 = 8.317e3; // N * m * kg^-1 * mol^-1 * K^-1
