/*
Interplanetary Planner
*/

use chrono::{DateTime, Utc};
use plotters::prelude::*;
use std::ops::Range;

use crate::donnager::{constants as cst, gravity::kepler as kepler, spacetime as xyzt};

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
    let mean_radius: f64 = (orb_radius_0 + orb_radius_f) / 2.;
    return (cst::SUN::GRAV_PARAM * (2. / orb_radius_0 - 1. / mean_radius)).sqrt();
}

/// Get solar system bodies
pub fn get_solar_system_bodies(
    planet_idxs: Vec<String>
) -> Vec<xyzt::Body> {
    let mut bodies: Vec<xyzt::Body>;

    // Mercury

    // Venus

    if planet_idxs.contains(&"earth".to_string()) {
        // Earth
        let earth: xyzt::Body = xyzt::Body {
            name: String::from("Earth"),
            grav_param: cst::EARTH::GRAV_PARAM,
            eq_radius: cst::EARTH::RADIUS_EQUATOR,
            rotation_rate: cst::EARTH::ROT_RATE,
            sidereal_day_hours: cst::EARTH::SIDEREAL_DAY,
            eccentricity: cst::EARTH::ECC
        };
        bodies.append(& mut vec![earth])
    }

    // Mars
    let mars: xyzt::Body = xyzt::Body{
        name: "Mars".to_string(),
        grav_param: cst::MARS::GRAV_PARAM,
        eq_radius: cst::MARS::RADIUS_EQUATOR,
        rotation_rate: cst::MARS::ROT_RATE,
        sidereal_day_hours: cst::MARS::SIDEREAL_DAY,
        eccentricity: cst::MARS::ECC
    };

    // Jupiter

    // Saturn

    // Uranus

    // Neptune
    let bodies = bodies;
    return bodies
}


/// Populate vector of solar system objects
pub fn get_solar_system_orbits(
    bodies: Vec<String>
) -> Vec<kepler::Orbit> {
    let mut orbits: Vec<kepler::Orbit>;

    // Sun
    let sun: xyzt::Body = xyzt::Body {
        name: "Sun".to_string(),
        grav_param: cst::SUN::GRAV_PARAM,
        eq_radius: cst::SUN::RADIUS_EQUATOR,
        rotation_rate: 0.,
        sidereal_day_hours: 0.,
        eccentricity: cst::SUN::ECC
    };

    let (year, month, day) = xyzt::julian_to_gregorian(
        cst::J2000_DAY as i32);
    let epoch_date_time = xyzt::ymd_hms_to_datetime(
        year, month as u32, day as u32, 0, 0, 0);

    // Mercury

    // Venus

    if bodies.contains(&"earth".to_string()) {
        // Earth-Sun orbit
        let earth_sun_orbit: kepler::Orbit = kepler::Orbit::from_keplerian(
            "Earth-Sun Orbit".to_string(),
            sun.clone(),
            cst::EarthSunOrbit::SEMI_MAJOR,
            cst::EarthSunOrbit::ECC,
            cst::EarthSunOrbit::INC,
            cst::EarthSunOrbit::RAAN,
            cst::EarthSunOrbit::ARG_PERIHELION,
            cst::EarthSunOrbit::MEAN_ANOMALY,
            cst::EarthSunOrbit::MEAN_MOTION,
            epoch_date_time);
        orbits.append(&mut vec![earth_sun_orbit]);
    }

    // Mars-Sun orbit
    let mars_sun_orbit: kepler::Orbit = kepler::Orbit::from_keplerian(
        "Mars-Sun Orbit".to_string(),
        sun.clone(),
        cst::MarsSunOrbit::SEMI_MAJOR,
        cst::MarsSunOrbit::ECC,
        cst::MarsSunOrbit::INC,
        cst::MarsSunOrbit::RAAN,
        cst::MarsSunOrbit::ARG_PERIHELION,
        cst::MarsSunOrbit::MEAN_ANOMALY,
        cst::MarsSunOrbit::MEAN_MOTION,
        epoch_date_time
    );

    
    // Jupiter

    // Saturn

    // Uranus

    // Neptune

    let orbits = orbits;
    return orbits
}

/// Calculate next hohmann transfer launch window
/// 
/// Inputs
/// ------
/// start_datetime: `DateTime<Utc>`
///     Start datetime for search
/// 
/// orbit_1: `Orbit`
///     Orbit of starting planet
/// 
/// orbit_2: `Orbit`
///     Orbit of ending planet
/// 
/// Outputs
/// -------
/// `Vec<DateTime<Utc>>`
///     Vector of launch windows
pub fn calc_next_hohmann_launch_window(
    start_datetime: DateTime<Utc>,
    orbit_1: kepler::Orbit,
    orbit_2: kepler::Orbit
) -> Vec<DateTime<Utc>>{

    let period_1 = orbit_1.calc_period();
    let period_2 = orbit_2.calc_period();
    let synodic_period = period_1 / period_2;

    let epoch_time = start_datetime.timestamp() as f64;
    let true_anonmaly_0_1 = orbit_1.calc_true_anomaly(epoch_time);

    let epoch_time = start_datetime.timestamp() as f64;
    let true_anonmaly_0_2 = orbit_2.calc_true_anomaly(epoch_time);

    let diff = true_anonmaly_0_2 - true_anonmaly_0_1;
    // Find when diff =  +/- 180
    let error = (180. - diff) / synodic_period;



}


pub fn show_trajectory(
    Orbits: Vec<kepler::Orbit>,
    Maneuvers: Vec<kepler::Maneuver>
){
    // Create drawing canvas of solar system bodies
    let pathname: String = format!("SolarSystemTrajectory.png");

    let drawing_area = 
        BitMapBackend::new(&pathname, (500, 500))
            .into_drawing_area();

    let mut chart_builder = ChartBuilder::on(&drawing_area);
    
    // Plot each orbital ellipse in the plane
    for orbit in Orbits{
        if orbit.central_body.name != "Sun"{
            // Calculate central body trajectory
            let planet_orbit = get_solar_system_orbits(
                vec![orbit.central_body.name.to_string()])[0];
            let motion_0 = planet_orbit.calc_motion(time, frame);
        } else {
            // Plot trajectory directly

        }
    }

    // Plot each maneuver location, magnitude, and direction


}

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
    start_date_time: DateTime<Utc>,
    stop_date_time: DateTime<Utc>,
    orbit_1: kepler::Orbit,
    orbit_2: kepler::Orbit
){
    let frame = xyzt::ReferenceFrames::InertialCartesian;
    
    for launch_time in start_date_time.timestamp()..stop_date_time.timestamp() {
        let motion1 = orbit_1.calc_motion(launch_time as f64, frame);
        let motion2 = orbit_2.calc_motion(launch_time as f64, frame);
        let v_inf: f64 = calc_esc_vel(motion1[0].norm(), motion2[0].norm());

    }

    //  calculate deltavs for two impulse maneuvers (min, max)
    //  calculate total time of flight for each

    let pathname: String = format!(
        "{}_to_{}_time_vs_fuel_{:?}.png",  
        orbit_1.central_body.name, 
        orbit_2.central_body.name,
        start_date_time);

    let plottitle: String = format!(
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
        (start_date_time.date_naive())..(stop_date_time.date_naive());  

    let longest_tof: chrono::Duration = 
        stop_date_time.date_naive() - start_date_time.date_naive();
    let y_spec: Range<chrono::Duration> = chrono::Duration{secs: 0, nanos: 0}..longest_tof;  

    let mut chart = 
            chart_builder.build_cartesian_2d(x_spec, y_spec).unwrap();
    chart
        .configure_mesh()
        .y_desc("Arrival Date")
        .x_desc("Departure Date")
        .draw()
        .unwrap();


    // Plot Delta v contours

    



    // Lines of -1 year slope
    for intercept in integer_tof_years:
        chart.draw_series(
            PointSeries::of_element(
                .iter().map(|p| (p.y, p.x)),
                1,
                &BLUE,
                &|c, s, st| {
                    Circle::new((c.0, c.1), s, st.filled())}
        )
    ).unwrap()
    .label("Const ToF")
    .legend(
        |(x, y)| 
        PathElement::new(vec![(x, y), (x + 20, y)], 
        &BLUE));
    

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

}
