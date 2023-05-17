/*
Aldrin Cycler between Earth and Mars
*/


use donnager::donnager::{constants as cst, spacetime as xyzt};


fn main() {
    // Configurable 
    // let start_date_time = xyzt::ymd_hms_to_datetime(2023, 1, 1, 12, 0, 0);
    // let stop_date_time = xyzt::ymd_hms_to_datetime(2053, 1, 1, 12, 0, 0);
   
    let start_alt: f64 = 408000.0;  
    let end_alt: f64 = 208000.0;  

    // let bodies = interplan::get_solar_system_bodies(
    //     vec!["earth".to_string(), "mars".to_string()]);

    // Givens
    // let planet_sun_orbits = interplan::get_solar_system_orbits(bodies);

    let semi_major_axis_earth: f64 = start_alt + cst::EARTH::RADIUS_EQUATOR;
    let semi_major_axis_mars: f64 = end_alt + cst::MARS::RADIUS_EQUATOR;

    // let launch_site: xyzt::SurfacePoint = xyzt::SurfacePoint {
    //     name: "Cape Canaveral Launch Site".to_string(),
    //     body: planets[0].clone(),
    //     pos_lla: Vector3::new(28.396837, -80.605659, 0.0)
    // };

    // let delta_v: f64 = launch_site.calc_delta_v(start_alt);

    // // Start in LEO
    // let orbit_one = grav::kepler::Orbit::from_keplerian(
    //     "Earth Parking".to_string(), 
    //     planets[0].clone(), 
    //     semi_major_axis_earth, 
    //     0., 
    //     0., 
    //     0., 
    //     0., 
    //     0., 
    //     grav::kepler::calc_mean_motion(
    //         semi_major_axis_earth, planets[0].grav_param), 
    //     start_date_time);

    // // End in LMO
    // let orbit_two = grav::kepler::Orbit::from_keplerian(
    //     "Mars Parking".to_string(), 
    //     planets[1].clone(), 
    //     semi_major_axis_mars, 
    //     0., 
    //     0., 
    //     0., 
    //     0., 
    //     0., 
    //     grav::kepler::calc_mean_motion(
    //         semi_major_axis_mars, planets[1].grav_param), 
    //         stop_date_time);
    
    
    // select start datetime

    // search forward for optimal launch windows

    // Compare patched conic vs 3bp fidelity

    // compare passive and active cyclers
    
    // plot trajectory, fuel, time

    // grav::interplan::show_porkchop_plots(
    //     start_date_time,
    //     stop_date_time,
    //     orbit_one,
    //     orbit_two
    // );

}
