/*
Aldrin Cycler between Earth and Mars
*/

use nalgebra::Vector3;
use chrono::DateTime;

use donnager::donnager::{constants as cst, spacetime as xyzt, gravity as grav};


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
        sun,
        cst::EarthSunOrbit::SEMI_MAJOR,
        cst::EarthSunOrbit::ORBIT_ECC,
        cst::EarthSunOrbit::INCLINATION,
        cst::EarthSunOrbit::LONG_ASC_NODE,
        cst::EarthSunOrbit::ARG_PERIAPSIS,
        cst::EarthSunOrbit::MEAN_ANOMALY
    );

    // Mars-Sun orbit
    let mars_orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_keplerian(
        "Mars_Sun",
        sun,
        cst::MarsSunOrbit::SEMI_MAJOR,
        cst::MarsSunOrbit::ECC,
        cst::MarsSunOrbit::INCLINATION,
        cst::MarsSunOrbit::LONG_ASC_NODE,
        cst::MarsSunOrbit::ARG_PERIAPSIS,
        cst::MarsSunOrbit::MEAN_ANOMALY
    );

 
    // select start datetime
    // search forward for optimal launch windows
    // Compare patched conic vs 3bp fidelity
    // compare passive and active cyclers
    // plot trajectory, fuel, time

    grav::interplan::show_porkchop_plots(
        (DateTime::parse_from_rfc3339("2020-01-01T00:00:00Z").unwrap(),
        DateTime::parse_from_rfc3339("2030-01-01T00:00:00Z").unwrap()),
        &earth_orbit,
        &mars_orbit
    )
}