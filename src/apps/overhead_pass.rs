use nalgebra::Matrix3;
use nalgebra::Vector3;
use sgp4;
use chrono::{DateTime, Utc};
use parse_tle::tle;

use donnager::donnager::spacetime as xyzt;
use donnager::donnager::constants as cst;

fn main() {

    // Parse CLI inputs
    // - Lat lon
    // - TLE of interest
    // - Search length

    let observer_lla: Vector3<f64> = Vector3::new(
        28.396837, 
        -80.605659, 
        0.0
    );

    // Earth
    let earth: xyzt::Body = xyzt::Body {
        name: String::from("Earth"),
        grav_param: cst::EARTH::GRAV_PARAM,
        eq_radius: cst::EARTH::RADIUS_EQUATOR,
        rotation_rate: cst::EARTH::ROT_RATE,
        sidereal_day_hours: cst::EARTH::SIDEREAL_DAY,
        eccentricity: cst::EARTH::SURFACE_ECC
    };
    
    let iss_str: &str = "
    ISS (ZARYA)
    1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
    2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";

    let now: DateTime<Utc> = Utc::now();
    let tle: tle::TLE = tle::parse(iss_str);

    let class_level: sgp4::Classification = match tle.classification.as_str() {
        "U" => sgp4::Classification::Unclassified,
        "C" => sgp4::Classification::Classified,
        "S" => sgp4::Classification::Secret,
        _ => sgp4::Classification::Unclassified
    };

    

    // NaiveDateTime = epoch_to_datetime(Epoch)



    let elements: sgp4::Elements = sgp4::Elements{
        object_name: Some(tle.name),
        international_designator: Some(tle.international_designator),
        norad_id: tle.catalog_number.parse::<u64>().expect("None numerical value found"),
        classification: class_level,
        datetime: tle.epoch,
        mean_motion_dot: tle.mean_motion_1,
        mean_motion_ddot: tle.mean_motion_2,
        drag_term: tle.radiation_pressure,
        element_set_number: tle.element_set_number,
        inclination: tle.inc,
        right_ascension: tle.raan,
        eccentricity: tle.eccentricity,
        argument_of_perigee: tle.arg_perigee,
        mean_anomaly: tle.mean_anomaly,
        mean_motion: tle.mean_motion,
        revolution_number: tle.rev_num as u64,
        ephemeris_type: tle.ephemeris_type,
    }; 

    let constants: sgp4::Constants = sgp4::Constants::from_elements(
        &elements
    ).unwrap();

    // Search next 24 hours?
    for min in 0..24*60 {
        println!("t = {} min", min);
        let minutes_since_epoch: sgp4::MinutesSinceEpoch = sgp4::MinutesSinceEpoch(min as f64);
        let prediction = constants.propagate(minutes_since_epoch).unwrap();
        let p_eci: [f64; 3] = prediction.position;
    }

    // TODO-TD: Vectorize and functionalize everything below for search space
    let p_eci: Vector3<f64> = Vector3::new(0., 0., 0.);

    let observer_ecef: Vector3<f64> = xyzt::planetodetic_to_cartesian_rotational(
        observer_lla, 
        cst::EARTH::RADIUS_EQUATOR, 
        cst::EARTH::SURFACE_ECC
    );

    let eci_to_ecef: Matrix3<f64> = xyzt::calc_inertial_rotational_rotam(
        now, 
        cst::EARTH::ROT_RATE * 60. * 60. * 24. 
    );

    let p_ecef: Vector3<f64> = eci_to_ecef * p_eci;
    let p_enu: Vector3<f64> = xyzt::fixed_frame_to_enu(
        observer_lla, 
        p_ecef, 
        cst::EARTH::RADIUS_EQUATOR, 
        cst::EARTH::SURFACE_ECC
    ); 

    let is_overhead: bool = p_enu[2] >= 0.;
    let observer_eci: Vector3<f64>  = eci_to_ecef.transpose() * observer_ecef;
    let is_night: bool = xyzt::is_eclipsed_by_earth(
        observer_eci, 
        now
    );
    let is_sunlit: bool = !xyzt::is_eclipsed_by_earth(
        p_eci, 
        now
    );
    let is_visible: bool = is_sunlit && is_overhead && is_night;
    if is_visible{
        let azelrad: Vector3<f64> = xyzt::enu_to_azelrad(p_enu);
        let nowstring: String = now.to_string();
        print!("{nowstring}");
    } else {
        // No sri
    }

}