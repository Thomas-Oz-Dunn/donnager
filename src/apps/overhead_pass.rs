use nalgebra::Vector3;
use sgp4;
use parse_tle::tle;

use donnager::donnager::spacetime as xyzt;
use donnager::donnager::constants as cst;

fn main() {

    let observer_lla: Vector3<f64> = Vector3::new(28.396837, -80.605659, 0.0);

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

    let tle_str = ;
    let tle: tle::TLE = tle::parse(tle_str);

    let elements: sgp4::Elements = 

    let constants: sgp4::Constants = sgp4::Constants::from_elements(&elements).unwrap();

    for hours in 0..24 {
        println!("t = {} min", hours * 60);
        let prediction = constants.propagate(sgp4::MinutesSinceEpoch((hours * 60) as f64)).unwrap();
        println!("    r = {:?} km", prediction.position);
        println!("    ṙ = {:?} km.s⁻¹", prediction.velocity);
    }


    let observer_ecef: Vector3<f64> = xyzt::planetodetic_to_cartesian_rotational(
        observer_lla, 
        cst::EARTH::RADIUS_EQUATOR, 
        cst::EARTH::SURFACE_ECC
    );

    let eci_to_ecef = xyzt::calc_inertial_rotational_rotam(
        date_time, 
        rot_rate_rad_day
    );
    let p_ecef: Vector3<f64> = eci_to_ecef * p_eci;
    let enu: Vector3<f64> = xyzt::fixed_frame_to_enu(
        observer_lla, 
        p_ecef, 
        cst::EARTH::RADIUS_EQUATOR, 
        cst::EARTH::SURFACE_ECC
    ); 

}