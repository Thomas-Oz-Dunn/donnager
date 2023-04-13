/*
Interplanetary Planner
*/

use nalgebra::{Vector3};
use chrono::{DateTime, Utc};
use plotters::prelude::*;
use std::ops::Range;

use crate::donnager::{constants as cst, gravity::kepler as kepler};

/// Calculate interplanetary mission escape velocity
/// 
/// Inputs
/// ------
/// orb_radius_0: `f64`
///     Starting planet heliocentric radius
/// 
/// orb_radius_f: `f64`
///     Ending planet heliocentric radius
pub fn calc_esc_vel(
    orb_radius_0: f64,
    orb_radius_f: f64,
) -> f64 {
    let av_radius: f64 = (orb_radius_0 + orb_radius_f) / 2.;
    return (cst::SUN::GRAV_PARAM *(2. / orb_radius_0 - 1. / av_radius)).sqrt();
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
//     // Find diff =  +/- 180
//     let distance = (180. - diff) / synodic_period;


// }


// pub fn show_trajectory(
//     Orbits: Vec<Orbit>,
//     Maneuvers: Vec<Maneuver>
// ){
//     // Create drawing canvas of solar system bodies

//     // Plot each orbital ellipse in the plane

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
pub fn show_porkchop_plots(
    datetime_lauch_window: (DateTime<Utc>, DateTime<Utc>),
    orbit_1: kepler::Orbit,
    orbit_2: kepler::Orbit
){

    // for launch_date in window
    //  calculate location of origin planet
    //  calculate deltavs for two impulse maneuvers (min, max)
    //  calculate total time of flight for each

    let pathname = format!(
        "{}_to_{}_time_vs_fuel_{:?}.png",  
        orbit_1.central_body.name, 
        orbit_2.central_body.name,
        datetime_lauch_window.0);

    let plottitle = format!(
        "{} to {}", 
        orbit_1.central_body.name, 
        orbit_2.central_body.name);
    let drawing_area = 
        BitMapBackend::new(&pathname, (500, 500))
            .into_drawing_area();

    let mut chart_builder = ChartBuilder::on(&drawing_area);

    drawing_area
        .fill(&WHITE)
        .unwrap();

    chart_builder
        .margin(5)
        .set_left_and_bottom_label_area_size(35)
        .caption(
            plottitle, (
                "Times New Roman", 
                20, 
                FontStyle::Bold, 
                &BLACK)
                .into_text_style(&drawing_area));
    let x_spec: Range<chrono::NaiveDate> = 
        (datetime_lauch_window.0.date_naive())..(datetime_lauch_window.1.date_naive());  

    // Determine from shortest trip from first date and longest trip from last date
    let y_spec: Range<f64> = shortest_tof..longest_tof;  

    let mut chart = 
            chart_builder.build_cartesian_2d(x_spec, y_spec).unwrap();
    chart
        .configure_mesh()
        .y_desc("Arrival Date")
        .x_desc("Departure Date")
        .draw()
        .unwrap();


    // Plot Delta v contours

    // Travel time lines

}


/// Calculate lagrange point locations in DU
/// 
/// Inputs
/// ------
/// mass_1: `f64`
///     Larger mass
/// 
/// mass_2: `f64`
///     Smaller mass
pub fn calc_lagrange_points(
    mass_1: f64,
    mass_2: f64
) -> Vec<Vector3<f64>> {
    let mass_ratio: f64= mass_2 / (mass_1 + mass_2);

    // Permissible if mass_2 is <<<< mass_1
    let xl12 = (mass_2/(3.*mass_1)).powf(1./3.);
    let xl3 = 7.*mass_2 / (12. * mass_1);

    let l1: Vector3<f64> = Vector3::new(1. - xl12, 0.,0.);
    let l2: Vector3<f64> = Vector3::new(1. + xl12 ,0.,0.);
    let l3: Vector3<f64> = Vector3::new(-xl3 ,0.,0.);
    let l4: Vector3<f64> = Vector3::new(mass_ratio - 0.5, -(3_f64).sqrt()/2., 0.);
    let l5: Vector3<f64> = Vector3::new(mass_ratio - 0.5, (3_f64).sqrt()/2., 0.);
    return vec![l1, l2, l3, l4, l5]
}


#[cfg(test)]
mod interplan_tests {

    use super::*;

    #[test]
    fn test_calc_esc_vel()
    {
        let orb_radius_0: f64 = 1.;
        let orb_radius_f: f64 = 2.;
        let esc_vel: f64 = calc_esc_vel(orb_radius_0, orb_radius_f);
        assert_eq!(esc_vel, 0.5);
    }

    #[test]
    fn test_calc_lagrange_points()
    {
        let mass_1: f64 = 1.;
        let mass_2: f64 = 0.001;
        let l_points = calc_lagrange_points(mass_1, mass_2);
        assert_eq!(l_points[0], Vector3::new(0.999, 0., 0.));
    }
}
