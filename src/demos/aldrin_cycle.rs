/*
Aldrin Cycler between Earth and Mars
*/

use nalgebra::Vector3;

use donnager::donnager::{constants as cst, spacetime as xyzt, gravity as grav};


fn main() {
    let (year, month, day) = xyzt::julian_to_gregorian(
        cst::J2000_DAY as i32);
    let epoch_date_time = xyzt::ymd_hms_to_datetime(
        year, month as u32, day as u32, 0, 0, 0);

    // Earth-Sun orbit
    let earth_orbit: grav::kepler::Orbit = grav::kepler::Orbit::from_keplerian(
        "Earth-Sun Orbit",
        xyzt::SUN.clone(),
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
        "Mars-Sun Orbit",
        xyzt::SUN.clone(),
        cst::MarsSunOrbit::SEMI_MAJOR,
        cst::MarsSunOrbit::ECC,
        cst::MarsSunOrbit::INC,
        cst::MarsSunOrbit::RAAN,
        cst::MarsSunOrbit::ARG_PERIHELION,
        cst::MarsSunOrbit::MEAN_ANOMALY,
        cst::MarsSunOrbit::MEAN_MOTION,
        epoch_date_time
    );

    
    // Start in LEO
    // TODO-TD: Calculate nominal LEO Orbit
    let pos_leo: Vector3<f64> = Vector3::new(x_pos_leo, y_pos_leo, z_pos_leo);
    let vel_leo: Vector3<f64> = Vector3::new(x_vel_leo, y_vel_leo, z_vel_leo);
    let orbit_one= grav::kepler::Orbit::from_pos_vel(
        "Earth Parking",
        xyzt::EARTH.clone(),
        pos_leo, 
        vel_leo,
        epoch_date_time
    );

    // End in LMO
    // TODO-TD: Calculate nominal LMO Orbit
    let pos_lmo: Vector3<f64> = Vector3::new(x_pos_lmo, y_pos_lmo, z_pos_lmo);
    let vel_lmo: Vector3<f64> = Vector3::new(x_vel_lmo, y_vel_lmo, z_vel_lmo);
    let orbit_two = grav::kepler::Orbit::from_pos_vel(
        "Mars Ending",
        xyzt::MARS.clone(),
        pos_lmo, 
        vel_lmos,
        epoch_date_time
    );
    // select start datetime
    // search forward for optimal launch windows
    // Compare patched conic vs 3bp fidelity
    // compare passive and active cyclers
    // plot trajectory, fuel, time

    let start_date_time = xyzt::ymd_hms_to_datetime(2023, 1, 1, 12, 0, 0);
    let stop_date_time = xyzt::ymd_hms_to_datetime(2053, 1, 1, 12, 0, 0);

    grav::interplan::show_porkchop_plots(
        start_date_time,
        stop_date_time,
        orbit_1,
        orbit_2
    );

}
