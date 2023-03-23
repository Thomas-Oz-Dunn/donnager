/*
Aldrin Cycler between Earth and Mars
*/

use nalgebra::Vector3;

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
        cst::EARTH_ORBIT_SEMI_MAJOR,
        cst::EARTH_ORBIT_ECC
    );

    // Mars-Sun orbit
    let mars_orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_keplerian(
        "Mars_Sun",
        sun::grav_param,
        cst::MARS_ORBIT_SEMI_MAJOR,
        cst::MARS_ORBIT_ECC
    );



    // select start datetime
    // search forward for optimal launch windows
    // compare passive and active cyclers
    // plot trajectory, fuel, time

}