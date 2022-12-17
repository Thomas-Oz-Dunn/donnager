/*
Constants

TODO-TD: query JPL SSDs
TODO-TD: load table of specific heats and molar weights
*/

// Space
pub const GRAV_CONST: f64 = 6.6743015e-11; // m^3 * kg^-1 * s^-2
pub const SUN_MASS: f64 = 1.9891e30; // kg

// Earth
pub const EARTH_MASS: f64 = 5.972e24; // kg
pub const EARTH_RADIUS_EQUATOR: f64 = 6.378137e6; // m
pub const EARTH_ECC: f64 = 0.08182;
pub const EARTH_ROT_RATE: f64 = 7.2921150e-5; // radians per second;
pub const EARTH_GRAV_PARAM: f64 = EARTH_MASS * GRAV_CONST; // m^3 * s^-2

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