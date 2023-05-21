/*
Constants
*/

// TODO-TD: Provide sources for every constant
use std::f64::consts::PI;

pub const CYCLES_TO_DEGREES: f64 = 360.;

/// Degrees to radians 
/// 
/// Source:
/// 
/// https://en.wikipedia.org/wiki/Radian
pub const DEG_TO_RAD: f64 = PI / 180.0;

/// Convert degrees to radians 
/// 
/// Source:
/// 
/// https://en.wikipedia.org/wiki/Radian
pub fn to_radians(degrees: f64) -> f64 {DEG_TO_RAD * degrees}

/// Radians to degrees 
/// 
/// Source:
/// 
/// https://en.wikipedia.org/wiki/Degree_(angle)
pub const RAD_TO_DEG: f64 = 180.0 / PI;

/// Convert radians to degrees
///   
/// Source:
/// 
/// https://en.wikipedia.org/wiki/Degree_(angle)
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


/// Newtonian Gravitational constant
/// 
/// Source:
/// https://en.wikipedia.org/wiki/Gravitational_constant
pub const GRAV_CONST: f64 = 6.6743015e-11; // m^3 * kg^-1 * s^-2

/// J2000 Reference Epoch Day
/// 
/// Source:
/// https://en.wikipedia.org/wiki/J2000
pub const J2000_DAY: f64 = 2451545.0;


/// Physical Constants of the Sun
/// 
/// Source:
/// https://en.wikipedia.org/wiki/Sun
pub struct SUN;
impl SUN {
    /// Mass of the Sun in kilograms
    pub const MASS: f64 = 1.9885e30; // kg
    /// Radius of the Sun in meters
    pub const RADIUS_EQUATOR: f64 = 6.957e8; // m
    /// Sun centered gravitational parameter
    pub const GRAV_PARAM: f64 = SUN::MASS * GRAV_CONST;
    pub const LUMINOSITY: f64 = 3.828e26; // Watts
    /// Mean Solar Flux at 1 AU
    pub const MEAN_SOLAR_FLUX: f64 = 1366.91; // Watts * m-2
    pub const ECC: f64 = 1.799991e-5;
}

/// Characteristics of Mercury
/// 
/// Source
/// ------
/// https://en.wikipedia.org/wiki/Mercury_(planet)
pub struct MERCURY;
impl MERCURY{
    pub const MASS: f64 = 3.302e23;
    pub const RADIUS_EQUATOR: f64 = 2439.7;
    pub const GRAV_PARAM: f64 = MERCURY::MASS * GRAV_CONST;
    pub const SURFACE_ECC: f64 = 7.004;
    pub const ORBIT_SEMI_MAJOR: f64 = 5.7909050e7; // km
    pub const ORBIT_ECC: f64 = 0.205630;
    pub const ARG_PERIHELION: f64 = 29.124 * DEG_TO_RAD;
    pub const INC: f64 = 7.005 * DEG_TO_RAD;
    pub const RAAN: f64 = 48.331 * DEG_TO_RAD;
    pub const MEAN_ANOMALY: f64 = 174.796 * DEG_TO_RAD; // radians
}

pub struct VENUS;
impl VENUS{
    pub const MASS: f64 = 4.869e24;
    pub const RADIUS_EQUATOR: f64 = 6051.8;
    pub const GRAV_PARAM: f64 = VENUS::MASS * GRAV_CONST;
    pub const ECC: f64 = 3.39471;
}

pub struct VenusSunOrbit;
impl VenusSunOrbit{
    // pub const SEMI_MAJOR: f64 = ; // km
}


/// Characteristics of Earth
/// 
/// Source
/// ------
/// https://en.wikipedia.org/wiki/Earth
pub struct EARTH;
impl EARTH {
    /// Mass of the Earth in kilograms
    pub const MASS: f64 = 5.972e24; // kg
    /// Radius of Earth at equator in meters
    pub const RADIUS_EQUATOR: f64 = 6.378137e6; // m
    pub const RADIUS_POLE: f64 = 6.3567e6;
    
    /// Earth centered gravitational parameter
    pub const GRAV_PARAM: f64 = EARTH::MASS * GRAV_CONST; // m^3 * s^-2
    pub const SURFACE_ECC: f64 = 0.08182;
    pub const SOLAR_DAY: f64 = 24.;
    pub const SIDEREAL_DAY: f64 = 23.9344696; // hours
    pub const ROT_RATE: f64 = 7.2921150e-5; // radians per second;
    pub const DAYS_PER_YEAR: f64 = 365.25;
    pub const AXIAL_TILT: f64 = -23.44;

    pub const SEMI_MAJOR: f64 = 149.60e6; // km
    pub const ORBIT_ECC: f64 = 0.0167086;
    pub const ARG_PERIHELION: f64 = 102.9372 * DEG_TO_RAD;
    pub const INC: f64 = 0.0 * DEG_TO_RAD;
    pub const RAAN: f64 = 0.0 * DEG_TO_RAD;
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


/// Characteristics of Mars
/// 
/// Source
/// ------
/// https://en.wikipedia.org/wiki/Mars

pub struct MARS;
impl MARS{
    pub const MASS: f64 = 6.4171e23; // kg
    pub const RADIUS_EQUATOR: f64 = 3.3895e6; // m
    pub const GRAV_PARAM: f64 = MARS::MASS * GRAV_CONST; // m^3 * s^-2
    pub const SURFACE_ECC: f64 = 0.0934;
    pub const SIDEREAL_DAY: f64 = 24.6229; // hours
    pub const ROT_RATE: f64 = 7.088e-5; // radians per second

    pub const ORBIT_SEMI_MAJOR: f64 = 2.27939366e8; // km
    pub const ORBIT_ECC: f64 = 0.0934;
    pub const ARG_PERIHELION: f64 = 286.5016 * DEG_TO_RAD;
    pub const INC: f64 = 1.850 * DEG_TO_RAD;
    pub const RAAN: f64 = 49.57854 * DEG_TO_RAD;
    pub const MEAN_MOTION: f64 = 0.524; // radians per day
    pub const MEAN_ANOMALY: f64 = 19.412 * DEG_TO_RAD; // radians
}

pub struct JUPITER;
impl JUPITER{
    pub const MASS: f64 = 1.899e27;
    pub const RADIUS_EQUATOR: f64 = 71492.; // m
    pub const GRAV_PARAM: f64 = MARS::MASS * GRAV_CONST; // m^3 * s^-2
    pub const ECC: f64 = 0.04838624;
    pub const SIDEREAL_DAY: f64 = 0.4135; // hours
    // pub const ROT_RATE: f64 = ; // radians per second
}

pub struct SATURN;
impl SATURN{
    pub const MASS: f64 = 5.685e26;
}

pub struct URANUS;
impl URANUS{
    pub const MASS: f64 = 8.682e25;
}

pub struct NEPTUNE;
impl NEPTUNE{
    pub const MASS: f64 = 1.024e26;
}


/// Speed of light in a vacuum `m * s^-1`
/// 
/// Source:
/// 
/// https://en.wikipedia.org/wiki/Speed_of_light
pub const SPEED_OF_LIGHT: f64 = 2.99792458e8;

pub const SOLAR_MASS_KG: f64 = 1.98892e30;

pub const DISTANCE_UNIT: f64 = 1.495978707e8; // meters

/// Electrical Constant F * m^-1
/// 
/// Source:
/// 
/// https://en.wikipedia.org/wiki/Vacuum_permittivity
pub const ELECTRIC_CONST: f64 = 8.8541878123e-12;

pub const WIENS_DIS_CONST: f64 = 2.897771955;

/// Planck's Constant J * Hz-1
/// 
/// Source:
/// 
/// https://en.wikipedia.org/wiki/Planck_constant
pub const PLANCKS_CONST: f64 = 6.62607015e-34;

/// Boltzman Constat J * K-1
/// 
/// Source:
/// 
/// https://en.wikipedia.org/wiki/Boltzmann_constant
pub const BOLTZMAN_CONST: f64 = 1.380649e-23; 

/// Mass of Proton kg
pub const MASS_PROTON: f64 = 1.67262e-27;
pub const MASS_NEUTRON: f64 = 1.67492749e-27; //kg
pub const ELECTRON_VOLT: f64 = 1.602e-19; // Joule


pub const MOLE: f64 = 6.02e23; // molecules / mole
pub const GAS_CONST: f64 = 8.317e3; // N * m * kg^-1 * mol^-1 * K^-1
