/*
Aldrin Cycler between Earth and Mars
*/

use nalgebra::Vector3;
use chrono::DateTime;

use donnager::donnager::constants as cst;
use donnager::donnager::spacetime as xyzt;
use donnager::donnager::gravity as grav;


fn main() {
    let earth: xyzt::Body = xyzt::Body {
        name: "Earth".to_string(),
        grav_param: cst::EARTH::GRAV_PARAM,
        eq_radius: cst::EARTH::RADIUS_EQUATOR,
        rotation_rate: cst::EARTH::ROT_RATE,
        eccentricity: cst::EARTH::ECC
    };

    // Mars
    let mars: xyzt::Body = xyzt::Body{
        name: "Mars".to_string(),
        grav_param: cst::MARS::GRAV_PARAM,
        eq_radius: cst::MARS::RADIUS_EQUATOR,
        rotation_rate: cst::MARS::ROT_RATE,
        eccentricity: cst::MARS::ECC
    };

    // Sun
    let sun: xyzt::Body = xyzt::Body {
        name: "Sun".to_string(),
        grav_param: cst::SUN::GRAV_PARAM,
        eq_radius: cst::SUN::RADIUS_EQUATOR,
        rotation_rate: cst::SUN::ROT_RATE,
        eccentricity: cst::SUN::ECC
    };

    // Earth-Sun orbit
    let earth_orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_keplerian(
        "Earth_Sun",
        sun::grav_param,
        cst::EarthSunOrbit::SEMI_MAJOR,
        cst::EarthSunOrbit::ORBIT_ECC
    );

    // Mars-Sun orbit
    let mars_orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_keplerian(
        "Mars_Sun",
        sun::grav_param,
        cst::MarsSunOrbit::SEMI_MAJOR,
        cst::MarsSunOrbit::ECC
    );
 
    // select start datetime
    let start_datetime = DateTime<UTC>

    // search forward for optimal launch windows
    // Compare patched conic vs 3bp fidelity
    // compare passive and active cyclers
    // plot trajectory, fuel, time

}