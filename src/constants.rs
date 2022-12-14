/*
Constants
*/

// Gravitational parameter
pub const GRAV_CONST: f64 = 6.6743015e-11; // m^3 * kg^-1 * s^-2

// Mass of the Sun in kilograms
pub const SUN_MASS: f64 = 1.9891e30; // kg

// Radius of the Sun in meters
pub const SUN_RADIUS_EQUATOR: f64 = 9.6e7; // m

// Sun centered gravitational parameter
pub const SUN_GRAV_PARAM: f64 = SUN_MASS * GRAV_CONST;

pub const SUN_LUMINOSITY: f64 = 3.828e26; // Watts

// Mass of the Earth in kilograms
pub const EARTH_MASS: f64 = 5.972e24; // kg

// Radius of Earth at equator in meters
pub const EARTH_RADIUS_EQUATOR: f64 = 6.378137e6; // m

// Earth centered gravitational parameter
pub const EARTH_GRAV_PARAM: f64 = EARTH_MASS * GRAV_CONST; // m^3 * s^-2

pub const EARTH_ECC: f64 = 0.08182;
pub const EARTH_ROT_RATE: f64 = 7.2921150e-5; // radians per second;

// Electromagnetism
pub const SPEED_OF_LIGHT: f64 = 2.99792458e8; // m * s^-1
pub const ELECTRIC_CONST: f64 = 8.8541878123; // F * m^-1
pub const PLANCKS_CONST: f64 = 6.62607015e-34; // J * Hz-1
pub const BOLTZMAN_CONST: f64 = 1.380649e-23; // J * K-1

// Nuclear physics
pub const MASS_PROTON: f64 = 1.67262e-27; //kg
pub const MASS_NEUTRON: f64 = 1.67492749e-27; //kg
pub const ELECTRON_VOLT: f64 = 1.602e-19; // Joule

// Chemistry
pub const MOLE: f64 = 6.02e23; // molecules / mole
pub const GAS_CONST: f64 = 8.317e3; // N * m * kg^-1 * mol^-1 * K^-1