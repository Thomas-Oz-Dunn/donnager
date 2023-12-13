use nalgebra::Vector;
use nalgebra::Vector3;
use sgp4;
use parse_tle::tle;

use donnager::donnager::spacetime as xyzt;
use donnager::donnager::constants as cst;

fn main() {

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


    // CLI inputs
    // - Datetime
    // - Lat lon
    // - TLE of interest
    let date_time: Datetime<Utc>;
    let tle_str: &str = "";
    let tle: tle::TLE = tle::parse(tle_str);
    let arg_peri: f64 = tle.arg_perigee;

    let elements: sgp4::Elements = sgp4::Elements{
        object_name: Some(tle.name),
        international_designator: Option<String>,
        norad_id: u64,
        classification: Classification,
        datetime: NaiveDateTime,
        mean_motion_dot: tle.mean_motion_1,
        mean_motion_ddot: tle.mean_motion_2,
        drag_term: tle.radiation_pressure,
        element_set_number: u64,
        inclination: tle.inc,
        right_ascension: tle.raan,
        eccentricity: tle.eccentricity,
        argument_of_perigee: tle.arg_perigee,
        mean_anomaly: tle.mean_anomaly,
        mean_motion: tle.mean_motion,
        revolution_number: u64,
        ephemeris_type: u8,
    }; 

    let constants: sgp4::Constants = sgp4::Constants::from_elements(
        &elements
    ).unwrap();

    let p_eci: Vector3<f64> = Vector3::new(0., 0., 0.);
    for min in 0..24*60 {
        println!("t = {} min", min);
        let minutes_since_epoch: sgp4::MinutesSinceEpoch = sgp4::MinutesSinceEpoch(min as f64);
        let prediction = constants.propagate(minutes_since_epoch).unwrap();
        let p_eci: [f64; 3] = prediction.position;
    }

    let is_bright: bool = xyzt::is_eclipsed_by_earth(
        p_eci, 
        date_time
    );
    let observer_ecef: Vector3<f64> = xyzt::planetodetic_to_cartesian_rotational(
        observer_lla, 
        cst::EARTH::RADIUS_EQUATOR, 
        cst::EARTH::SURFACE_ECC
    );

    let eci_to_ecef = xyzt::calc_inertial_rotational_rotam(
        date_time, 
        cst::EARTH::ROT_RATE * 60. * 60. * 24. 
    );

    let p_ecef: Vector3<f64> = eci_to_ecef * p_eci;
    let enu: Vector3<f64> = xyzt::fixed_frame_to_enu(
        observer_lla, 
        p_ecef, 
        cst::EARTH::RADIUS_EQUATOR, 
        cst::EARTH::SURFACE_ECC
    ); 

}