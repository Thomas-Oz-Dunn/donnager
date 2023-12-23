use sgp4;
use nalgebra::{Matrix3, Vector3};
use chrono::{DateTime, Utc, NaiveDate, NaiveTime, NaiveDateTime, Duration};
use parse_tle::tle;
use clap::{Parser, Args};

use donnager::donnager::{xyzt, constants as cst};

pub fn calc_overhead_passes(
    observer_lla: Vector3<f64>,
    days_to_search: i64,
    object_name: Option<String>,
    international_designator: Option<String>,
    epoch: NaiveDateTime,
    mean_motion_dot: f64,
    mean_motion_ddot: f64,
    drag_term: f64,
    element_set_number: u64,
    inclination: f64,
    right_ascension: f64,
    eccentricity: f64,
    argument_of_perigee: f64,
    mean_anomaly: f64,
    mean_motion: f64,
    revolution_number: u64,
    ephemeris_type: u8,
) -> (Vec<Vector3<f64>>, Vec<DateTime<Utc>>) {

    let elements: sgp4::Elements = sgp4::Elements{
        object_name: object_name,
        international_designator: international_designator,
        norad_id: 0,
        classification: sgp4::Classification::Unclassified,
        datetime: epoch,
        mean_motion_dot: mean_motion_dot,
        mean_motion_ddot: mean_motion_ddot,
        drag_term: drag_term,
        element_set_number: element_set_number,
        inclination: inclination,
        right_ascension: right_ascension,
        eccentricity: eccentricity,
        argument_of_perigee: argument_of_perigee,
        mean_anomaly: mean_anomaly,
        mean_motion: mean_motion,
        revolution_number: revolution_number,
        ephemeris_type: ephemeris_type,
    }; 

    let constants: sgp4::Constants = sgp4::Constants::from_elements(
        &elements
    ).unwrap();

    let now: DateTime<Utc> = Utc::now();
    
    let observer_ecef: Vector3<f64> = xyzt::planetodetic_to_cartesian_rotational(
        observer_lla, 
        cst::EARTH::RADIUS_EQUATOR, 
        cst::EARTH::SURFACE_ECC
    );
    
    let min_observer: i64 = (now - epoch.and_utc()).num_minutes();
    let min_dur: i64 = days_to_search * 24 * 60;
    
    // ^ setup

    // v run

    // TODO-TD: Use .map()
    let mut azelrad: Vec<Vector3<f64>> = vec![];
    let mut times: Vec<DateTime<Utc>> = vec![];

    let now: DateTime<Utc> = Utc::now();
    
    let observer_ecef: Vector3<f64> = xyzt::planetodetic_to_cartesian_rotational(
        observer_lla, 
        cst::EARTH::RADIUS_EQUATOR, 
        cst::EARTH::SURFACE_ECC
    );

    for minut in 0..min_dur {
        let eval_date_time: DateTime<Utc> = now + Duration::minutes(minut);
        let min_since_epoch: sgp4::MinutesSinceEpoch = sgp4::MinutesSinceEpoch((minut + min_observer) as f64);
        let pos: [f64; 3] = constants.propagate(min_since_epoch).unwrap().position;
        let p_eci: Vector3<f64> = Vector3::<f64>::new(pos[0], pos[1], pos[2]);

        let eci_to_ecef: Matrix3<f64> = xyzt::calc_inertial_rotational_rotam(
            eval_date_time, 
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
            eval_date_time
        );

        let is_sunlit: bool = !xyzt::is_eclipsed_by_earth(
            p_eci, 
            eval_date_time
        );

        if is_sunlit && is_overhead && is_night{
            let azelra: Vector3<f64> = xyzt::enu_to_azelrad(p_enu);
            azelrad.append(&mut vec![azelra]);
            times.append(&mut vec![eval_date_time]);
        }
    }

    let results = (azelrad, times);
    return results;
}


fn main() {

    // Parse CLI inputs
    // - Lat lon
    // - TLE of interest
    // - Search length

    let observer_lla: Vector3<f64> = Vector3::new(
        28.396837, 
        -80.605659, 
        10.
    );
    
    let sample_tle: &str = 
    "POLAR                   
    1 23802U 96013A   23343.66857320  .00000239  00000+0  00000+0 0  9998
    2 23802  78.9918 236.4491 5800292 240.2923  50.2600  1.29847017132903";

    let days_to_search: i64 = 10;
    let is_verbose: bool = true;

    // ^ parse cli

    // v setup
    let tle: tle::TLE = tle::parse(sample_tle);

    let class_level: sgp4::Classification = match tle.classification.as_str() {
        "U" => sgp4::Classification::Unclassified,
        "C" => sgp4::Classification::Classified,
        "S" => sgp4::Classification::Secret,
        _ => sgp4::Classification::Unclassified
    };

    let (y, m, d, h, mi, s, milli) = tle.epoch.to_gregorian_utc();
    
    let nd: NaiveDate = NaiveDate::from_ymd_opt(
        y, 
        m as u32, 
        d as u32
    ).unwrap();

    let t: NaiveTime = NaiveTime::from_hms_milli_opt(
        h as u32, 
        mi as u32, 
        s as u32, 
        milli
    ).unwrap();

    let epoch: NaiveDateTime = NaiveDateTime::new(nd, t);       
    let passes = calc_overhead_passes(
        observer_lla,
        days_to_search,
        Some(tle.name),
        Some(tle.international_designator),
        epoch,
        tle.mean_motion_1,
        tle.mean_motion_2,
        tle.radiation_pressure,
        tle.element_set_number,
        tle.inc,
        tle.raan,
        tle.eccentricity,
        tle.arg_perigee,
        tle.mean_anomaly,
        tle.mean_motion,
        tle.rev_num as u64,
        tle.ephemeris_type,
    );

    if is_verbose{
        for (pos, datetime) in passes.0.iter().zip(passes.1.iter()){
            print!(datetime.as_str());
        }
    };

}
