/*
Aldrin Cycler between Earth and Mars
*/

use nalgebra::Vector3;

use donnager::donnager::{constants as cst, spacetime as xyzt, gravity::{self as grav, interplan}};


fn main() {
    // Configurable 
    let start_date_time = xyzt::ymd_hms_to_datetime(2023, 1, 1, 12, 0, 0);
    let stop_date_time = xyzt::ymd_hms_to_datetime(2053, 1, 1, 12, 0, 0);
   
    // TODO-TD: Calculate nominal LEO Orbit
    let earth_alt: f64 = 408000.0;  // LEO
    let pos_earth: Vector3<f64> = Vector3::new(x_pos_leo, y_pos_leo, z_pos_leo);
    let vel_earth: Vector3<f64> = Vector3::new(x_vel_leo, y_vel_leo, z_vel_leo);

    // TODO-TD: Calculate nominal LMO Orbit
    let pos_mars: Vector3<f64> = Vector3::new(x_pos_lmo, y_pos_lmo, z_pos_lmo);
    let vel_mars: Vector3<f64> = Vector3::new(x_vel_lmo, y_vel_lmo, z_vel_lmo);


    // Givens
    let planets = interplan::get_solar_system_bodies();
    let planet_orbits = interplan::get_solar_system_orbits(start_date_time);
    
    // Start in LEO
    let orbit_one= grav::kepler::Orbit::from_pos_vel(
        "Earth Parking".to_string(),
        EARTH.clone(),
        pos_earth, 
        vel_earth,
        epoch_date_time
    );

    // End in LMO
    let orbit_two = grav::kepler::Orbit::from_pos_vel(
        "Mars Ending".to_string(),
        MARS.clone(),
        pos_mars, 
        vel_mars,
        epoch_date_time
    );
    
    let delta_v: f64 = launch_site.calc_delta_v(altitude);
    let dev = orbit_two.
    // select start datetime
    // search forward for optimal launch windows
    // Compare patched conic vs 3bp fidelity
    // compare passive and active cyclers
    // plot trajectory, fuel, time

    grav::interplan::show_porkchop_plots(
        start_date_time,
        stop_date_time,
        orbit_1,
        orbit_2
    );

}
