/*
Aldrin Cycler between Earth and Mars
*/

use donnager::donnager::{constants as cst, spacetime as xyzt, gravity as grav};


fn main() {
    // Earth
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
        rotation_rate: 0.,
        eccentricity: cst::SUN::ECC
    };

    let (year, month, day) = xyzt::julian_to_gregorian(
        cst::J2000_DAY as i32);
    let epoch_date_time = xyzt::ymd_hms_to_datetime(
        year, month as u32, day as u32, 0, 0, 0);

    // Earth-Sun orbit
    let earth_orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_keplerian(
        "Earth-Sun Orbit".to_string(),
        sun.clone(),
        cst::EarthSunOrbit::SEMI_MAJOR,
        cst::EarthSunOrbit::ECC,
        cst::EarthSunOrbit::INC,
        cst::EarthSunOrbit::RAAN,
        cst::EarthSunOrbit::ARG_PERIHELION,
        cst::EarthSunOrbit::MEAN_ANOMALY,
        cst::EarthSunOrbit::MEAN_MOTION,
        epoch_date_time
    );

    // Mars-Sun orbit
    let mars_orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_keplerian(
        "Mars-Sun Orbit".to_string(),
        sun,
        cst::MarsSunOrbit::SEMI_MAJOR,
        cst::MarsSunOrbit::ECC,
        cst::MarsSunOrbit::INC,
        cst::MarsSunOrbit::RAAN,
        cst::MarsSunOrbit::ARG_PERIHELION,
        cst::MarsSunOrbit::MEAN_ANOMALY,
        cst::MarsSunOrbit::MEAN_MOTION,
        epoch_date_time
    );

 
    // select start datetime
    // search forward for optimal launch windows
    // Compare patched conic vs 3bp fidelity
    // compare passive and active cyclers
    // plot trajectory, fuel, time

    let start_date_time = xyzt::ymd_hms_to_datetime(2023, 1, 1, 12, 0, 0);
    let stop_date_time = xyzt::ymd_hms_to_datetime(2053, 1, 1, 12, 0, 0);
    let window = (start_date_time, stop_date_time);

    grav::interplan::show_porkchop_plots(
        window,
        earth_orbit,
        mars_orbit
    )

}