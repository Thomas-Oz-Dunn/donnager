/*
Constants

TODO-TD: query JPL SSDs
TODO-TD: load table of specific heats and molar weights
*/

// Space
pub const SUN_MASS: f64 = 1.9891e30; // kg
pub const EARTH_MASS: f64 = 5.972e24; // kg
pub const EARTH_RADIUS_EQUATOR: f64 = 6.378137e6; // m
pub const EARTH_ECC: f64 = 0.08182;

// Physics
pub const GRAV_CONST: f64 = 6.6743015e-11; // m^3 * kg^-1 * s^-2
pub const GAS_CONST: f64 = 8.317e3; // N * m * kg^-1 * mol^-1 * K^-1
pub const SPEED_OF_LIGHT: f64 = 2.99792458e8; // m * s^-1
pub const MOLE: f64 = 6.02e23; // molecules / mole