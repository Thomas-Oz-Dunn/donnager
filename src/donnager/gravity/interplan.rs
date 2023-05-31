/*
Interplanetary Mission Planner
*/
use polars::prelude::*;
use chrono::{DateTime, Utc};
use plotters::prelude::*;
use std::ops::Range;

use crate::donnager::{
    constants as cst, 
    gravity::kepler as kepler, 
    spacetime as xyzt};

/// Calculate interplanestary mission escape velocity
/// 
/// Inputs
/// ------
/// orb_radius_0: `f64`
///     Starting planet heliocentric radius
/// 
/// orb_radius_f: `f64`
///     Ending planet heliocentric radius
pub fn calc_esc_vel(
    grav_param: f64,
    orb_radius_0: f64,
    orb_radius_f: f64,
) -> f64 {
    let mean_radius: f64 = (orb_radius_0 + orb_radius_f) / 2.;
    return (grav_param * (2. / orb_radius_0 - 1. / mean_radius)).sqrt();
}

/// Get ephemeris dataframe
pub fn get_ephemeris() -> DataFrame {

    let ephemeris = df! (
        "number" => &[1, 2, 3, 4, 5, 6, 7, 8],
        "names" => &["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"],
        "groups" => &["Inner", "Inner", "Inner", "Inner", "Outer", "Outer", "Outer", "Outer"],
        "type" => &["Rocky", "Rocky", "Rocky", "Rocky", "Gas Giant", "Gas Giant", "Ice Giant", "Ice Giant"],
        "mass" => &[
            cst::MERCURY::MASS, 
            cst::VENUS::MASS, 
            cst::EARTH::MASS, 
            cst::MARS::MASS, 
            cst::JUPITER::MASS, 
            cst::SATURN::MASS, 
            cst::URANUS::MASS, 
            cst::NEPTUNE::MASS],

        "orbit_semi_major" => &[
            cst::MercurySunOrbit::SEMI_MAJOR, 
            cst::VenusSunOrbit::SEMI_MAJOR, 
            cst::EarthSunOrbit::SEMI_MAJOR, 
            cst::MarsSunOrbit::SEMI_MAJOR, 
            None, 
            None, 
            None, 
            None],

        "orbit_ecc" => &[
            cst::MercurySunOrbit::ECC, 
            cst::VenusSunOrbit::ECC, 
            cst::EarthSunOrbit::ECC, 
            cst::MarsSunOrbit::ECC, 
            None, 
            None, 
            None, 
            None],
        
        "argument_perigelion" => &[
            cst::MercurySunOrbit::ARG_PERIHELION, 
            cst::VenusSunOrbit::ARG_PERIHELION, 
            cst::EarthSunOrbit::ARG_PERIHELION, 
            cst::MarsSunOrbit::ARG_PERIHELION, 
            None, 
            None, 
            None, 
            None],

    ).unwrap();

    return ephemeris
}



// /// Calculate next hohmann transfer launch window
// /// 
// /// Inputs
// /// ------
// /// start_datetime: `DateTime<Utc>`
// ///     Start datetime for search
// /// 
// /// orbit_1: `Orbit`
// ///     Orbit of starting planet
// /// 
// /// orbit_2: `Orbit`
// ///     Orbit of ending planet
// /// 
// /// Outputs
// /// -------
// /// `Vec<DateTime<Utc>>`
// ///     Vector of launch windows
// pub fn calc_next_hohmann_launch_window(
//     start_datetime: DateTime<Utc>,
//     orbit_1: kepler::Orbit,
//     orbit_2: kepler::Orbit
// ) -> Vec<DateTime<Utc>>{

//     let period_1 = orbit_1.calc_period();
//     let period_2 = orbit_2.calc_period();
//     let synodic_period = period_1 / period_2;

//     let epoch_time = start_datetime.timestamp() as f64;
//     let true_anonmaly_0_1 = orbit_1.calc_true_anomaly(epoch_time);

//     let epoch_time = start_datetime.timestamp() as f64;
//     let true_anonmaly_0_2 = orbit_2.calc_true_anomaly(epoch_time);

//     let diff = true_anonmaly_0_2 - true_anonmaly_0_1;
//     // Find when diff =  +/- 180
//     let error = (180. - diff) / synodic_period;



// }


// pub fn show_trajectory(
//     Orbits: Vec<kepler::Orbit>,
//     Maneuvers: Vec<kepler::Maneuver>
// ){
//     // Create drawing canvas of solar system bodies
//     let pathname: String = format!("SolarSystemTrajectory.png");

//     let drawing_area = 
//         BitMapBackend::new(&pathname, (500, 500))
//             .into_drawing_area();

//     let mut chart_builder = ChartBuilder::on(&drawing_area);
    
//     // Plot each orbital ellipse in the plane
//     for orbit in Orbits{
//         if orbit.central_body.name != "Sun"{
//             // Calculate central body trajectory
//             let planet_orbit = get_solar_system_orbits(
//                 vec![orbit.central_body.name.to_string()])[0];
//             let motion_0 = planet_orbit.calc_motion(time, frame);
//         } else {
//             // Plot trajectory directly

//         }
//     }

//     // Plot each maneuver location, magnitude, and direction


// }

/// Show orbital transfer porkchop plots
/// 
/// Inputs
/// ------
/// datetime_launch_window: `(DateTime<Utc>, DateTime<Utc>)`
///     Start and end datetime for plot
/// 
/// orbit_1: `Orbit`
///     Orbit of starting planet
/// 
/// orbit_2: `Orbit`
///     Orbit of ending planet
pub fn calc_porkchop_plots(
    start_date_time: DateTime<Utc>,
    stop_date_time: DateTime<Utc>,
    orbit_1: kepler::Orbit,
    orbit_2: kepler::Orbit
){
    let frame = xyzt::ReferenceFrames::InertialCartesian;
    
    for launch_time in start_date_time.timestamp()..stop_date_time.timestamp() {
        let motion1 = orbit_1.calc_motion(launch_time as f64, frame);
        let motion2 = orbit_2.calc_motion(launch_time as f64, frame);
        let v_inf: f64 = calc_esc_vel(
            cst::SUN::GRAV_PARAM,
            motion1[0].norm(), 
            motion2[0].norm());

        }
    let lambert = lambert(orb_dpt, orb_arr);
            
    // Get norm delta velocities
    let dv_dpt = nalgebra.norm(man_lambert.impulses[0][1]);
    let dv_arr = nalgebra.norm(man_lambert.impulses[1][1]);
            
    // Compute all the output variables
    let c3_launch = dv_dpt**2;
    let c3_arrival = dv_arr**2;

    


}


#[cfg(test)]
mod interplan_tests {

    use super::*;

    #[test]
    fn test_calc_esc_vel()
    {
        let grav_param: f64 = cst::SUN::GRAV_PARAM;
        let orb_radius_0: f64 = 1000000.;
        let orb_radius_f: f64 = 2000000.;
        let esc_vel: f64 = calc_esc_vel(grav_param, orb_radius_0, orb_radius_f);
        assert_eq!(esc_vel, 13302555.410020536);
    }


}
