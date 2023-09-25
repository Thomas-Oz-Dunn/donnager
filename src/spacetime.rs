/*
Space time related functions
 */

use nalgebra::{Vector3, Matrix3};
use std::f64::consts::PI;
use chrono::prelude::*;
// use polars::prelude::*;

use crate::constants as cst;

/// Gravitational Body
/// 
/// Properties
/// ----------
/// name: `String`
///     Name of body
/// 
/// grav_param: `f64`
///     Gravitational Parameter
/// 
/// eq_radius: `f64`
///     Equatorial Radius in km
/// 
/// rotation_rate: `f64`
///     Rotation rate  in rad /s
/// 
/// sidereal_day_hours: `f64`
///     Length of sidereal day
/// 
/// eccentricity: `f64`
///     Polar Equatorial Eccentricity
#[derive(Clone, Debug, PartialEq)]
pub struct Body{
    pub name: String,
    pub grav_param: f64,
    pub eq_radius: f64,
    pub rotation_rate: f64,
    pub sidereal_day_hours: f64,
    pub eccentricity: f64
}


/// Reference Frames
/// 
/// Valid Frames
/// ------------
/// InertialCartesian 
///     Ex: (ECI)
/// RotationalCartesian 
///     Ex: (ECEF)
/// Perifocal 
///     Ex: (PFCL)
/// Planetodetic 
///     Ex: (Geodetic)
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ReferenceFrames {
    InertialCartesian, 
    RotationalCartesian,
    Perifocal,
    Planetodetic
}

/// Gravitational Body
impl Body {

    /// Convert into simple particle type
    /// 
    /// Inputs
    /// ------
    /// motion: `Vec<Vector3<f64>>`
    ///     motion of body in frame
    /// 
    /// Outputs
    /// -------
    /// particle: `Particle`
    ///     `Particle` equivalent of body
    pub fn to_particle(&self, motion: Vec<Vector3<f64>>) -> Particle {
        let mass: f64 = self.grav_param / cst::GRAV_CONST;
        let particle: Particle = Particle {mass: mass.clone(), motion: motion};
        return particle
    }


    /// Calculate gravitational acceleration at radial distance from Body
    /// 
    /// Inputs
    /// ------
    /// radius: `Vector3<f64>`
    ///     Radius in km from Body center
    /// 
    /// Outputs
    /// -------
    /// grav_acc: `Vector3<f64>`
    ///     Acceleration rate due to gravity
    pub fn calc_body_grav(&self, radius: Vector3<f64>) -> Vector3<f64> {
        let grav_acc: Vector3<f64>  = self.grav_param * radius / radius.norm().powi(3);
        return grav_acc
    }
    

    /// Calculate required orbital velocity at radial distance
    /// 
    /// Inputs
    /// ------
    /// radius: `f64`
    ///     Radius in km from Body center
    ///     
    /// Outputs
    /// -------
    /// val: `f64`
    ///     Required tangential velocity magnitude
    pub fn calc_orbital_velocity_mag(&self, radius: f64) -> f64 {
        let vel: f64 = (self.grav_param / radius).sqrt();
        return vel
    }


    /// Calculate required escape velocity at radial distance
    /// 
    /// Inputs
    /// ------
    /// radius: `f64`
    ///     Radius in km from Body center
    ///     
    /// Outputs
    /// -------
    /// val: `f64`
    ///     Required tangential velocity magnitude
    pub fn calc_escape_velocity_mag(&self, radius: f64) -> f64 {
        let vel: f64 = (2f64).sqrt() * self.calc_orbital_velocity_mag(radius);
        return vel
    }


    /// Calculate radius for stationary orbit above body surface
    /// 
    /// Outputs
    /// -------
    /// r_mag: `f64`
    ///     Magnitude of radius for stationary orbit
    pub fn calc_stationary_orbit(&self) -> f64 {
        let period: f64 = 2.0 * PI / self.rotation_rate;
        let semi_major: f64 = self.grav_param * period.powi(2); 
        let r_mag: f64 = (semi_major / (4.0 * PI.powi(2))).powf(1.0 / 3.0); 
        return r_mag
    }

}

/// Gravitational Particle
#[derive(Clone, Debug, PartialEq)]
pub struct Particle {
    pub mass: f64,
    pub motion: Vec<Vector3<f64>>,
}

impl Particle {

    /// Convert into Body type
    /// 
    /// Inputs
    /// ------
    /// name : `String`
    ///     Name of body
    /// 
    /// eq_radius : `f64`
    ///     Equatorial radius of body
    /// 
    /// rotation_rate : `f64`
    ///     Rotation rate of body
    /// 
    /// eccentricity : `f64`
    ///     Body oblateness
    pub fn to_body(
        &self, 
        name: String, 
        eq_radius: f64, 
        rotation_rate: f64, 
        sidereal_day_hours: f64,
        eccentricity: f64
    ) -> Body {
        let grav_param: f64 = self.mass * cst::GRAV_CONST;
        let body: Body = Body {
            name,
            grav_param,
            eq_radius,
            rotation_rate,
            sidereal_day_hours,
            eccentricity
        };
        body
    }
}


#[derive(Clone, Debug, PartialEq)]
pub struct SurfacePoint{
    pub name: String,
    pub body: Body,
    pub pos_lla: Vector3<f64>,
}


/// SurfacePoint methods
impl SurfacePoint{
    /// Calculate tangential velocity magnitude at surface point
    /// 
    /// Outputs
    /// -------
    /// tan_vel : `f64`
    ///     Tangential velocity magnitude in km / s
    pub fn calc_surface_vel(&self) -> f64 {
        let radius_km: f64 = self.body.eq_radius / cst::KILO;
        let equatorial_vel: f64 = self.body.rotation_rate * radius_km;
        let tan_vel: f64 = (self.pos_lla[0].cos() * equatorial_vel).abs();
        return tan_vel
    }
    
    /// Calculate radius from center of body at surface point
    /// 
    /// Outputs
    /// -------
    /// radius : `Vector3<f64>`
    ///     Radial vector in meters
    pub fn calc_surface_radius(&self) -> Vector3<f64> {
        let prime_vertical: f64 = calc_prime_vertical(
            self.pos_lla[0], 
            self.body.eq_radius,
            self.body.eccentricity);
        let dir: Vector3<f64> = planetodetic_to_cartesian_rotational(
            self.pos_lla,
            self.body.eq_radius,
        self.body.eccentricity);
        let radius: Vector3<f64> = (prime_vertical + self.pos_lla[2]) * dir;
        return radius
    }

    /// Calculate delta v required to reach an altitude from surface
    /// 
    /// Inputs
    /// ------
    /// altitude : `f64`
    ///     Altitude off of surface
    pub fn calc_delta_v(
        &self,
        altitude: f64,
    ) -> f64 {
        let radius: f64 = altitude + self.calc_surface_radius().norm();
        let surface_vel: f64 = self.calc_surface_vel();
        let delta_v_req: f64 = self.body.calc_orbital_velocity_mag(radius); 
        let delta_v: f64 = delta_v_req - surface_vel;
        return delta_v
    }


}

// /// Get ephemeris dataframe
// pub fn get_ephemeris() -> DataFrame {

//     let ephemeris = df! (
//         "number" => &[1, 2, 3, 4, 5, 6, 7, 8],
//         "name" => &[
//             "Mercury", "Venus", "Earth", "Mars", 
//             "Jupiter", "Saturn", "Uranus", "Neptune"],
//         "group" => &[
//             "Inner", "Inner", "Inner", "Inner", 
//             "Outer", "Outer", "Outer", "Outer"],
//         "type" => &[
//             "Rocky", "Rocky", "Rocky", "Rocky", 
//             "Gas Giant", "Gas Giant", "Ice Giant", "Ice Giant"],
//         "mass" => &[
//             cst::MERCURY::MASS, 
//             cst::VENUS::MASS, 
//             cst::EARTH::MASS, 
//             cst::MARS::MASS, 
//             cst::JUPITER::MASS, 
//             cst::SATURN::MASS, 
//             cst::URANUS::MASS, 
//             cst::NEPTUNE::MASS
//         ],
//         "equatorial_radius" => &[
//             cst::MERCURY::RADIUS_EQUATOR,
//             cst::VENUS::RADIUS_EQUATOR,
//             cst::EARTH::RADIUS_EQUATOR,
//             cst::MARS::RADIUS_EQUATOR,
//             cst::JUPITER::RADIUS_EQUATOR,
//             cst::SATURN::RADIUS_EQUATOR,
//             cst::URANUS::RADIUS_EQUATOR,
//             cst::NEPTUNE::RADIUS_EQUATOR
//         ],
//     ).unwrap();

//     return ephemeris
// }


/// Calculate inertial to corotational frame rotation matrix
/// 
/// Inputs
/// ------
/// date_time: `DateTime<Utc>`
///     Date and Time 
/// 
/// rot_rate_rad_day: `f64`
///     Rotation rate of body in radians per day
pub fn calc_inertial_rotational_rotam(
    date_time: DateTime<Utc>,
    rot_rate_rad_day: f64
) -> Matrix3<f64> {
    let year = date_time.year();
    let month = date_time.month() as i32;
    let day = date_time.day() as i32;

    let hours = date_time.hour() as f64;
    let minutes = date_time.minute() as f64;
    let seconds = date_time.second() as f64;

    let julian_day: f64 = date_to_julian_day_num(year, month, day) as f64;

    let sidereal_time: f64 = (hours + (minutes + seconds / 60.) / 60.) / 24.;
    let theta: f64 = 
        rot_rate_rad_day * (julian_day + sidereal_time - cst::J2000_DAY);

    let rotam: Matrix3<f64> = Matrix3::<f64>::new(
        theta.cos(), -theta.sin(), 0.,
        theta.sin(), theta.cos(), 0.,
        0., 0., 1.
    );

    return rotam
}

/// Calculate hours of sunlight at lattitude
/// 
/// Inputs
/// ------
/// lat_deg : `f64`
///     Latitude of coordinates in degrees
/// 
/// long_deg : `f64`
///     Longitude of coordinates in degrees
/// 
/// julian_day : `i32`
///     julian_date
pub fn calc_earth_day_length(
    lattitude_deg: f64,
    longitude_deg: f64,
    julian_day: i32
) -> f64 {
    let cycle_to_deg = cst::CYCLES_TO_DEGREES;
    let cycle_to_rad = cycle_to_deg * cst::RAD_TO_DEG;
    let ax_tilt = cst::EARTH::AXIAL_TILT;

    let j2000_days: f64 = (julian_day as f64) - cst::J2000_DAY;
    let j_star: f64 = j2000_days - longitude_deg / cycle_to_deg;
    
    // FIXME-TD: remove `magic numbers`
    let m: f64 = (357.5291 + 0.98560028 * j_star) % cycle_to_deg;
    let equat_center: f64 = 1.9148*(
        m.sin()) + 0.02*((2.*m).sin()) + 0.0003*((3.*m).sin()
    );
    let earth_tilt_adjust: f64 = (
        m + equat_center + cycle_to_deg / 2. + 
        cst::EARTH::ARG_PERIHELION) % cycle_to_deg;

    let sun_dec_sin: f64 = earth_tilt_adjust.sin() * ax_tilt.sin();
    let sun_dec_rad: f64 = sun_dec_sin.asin();

    let tan_lat: f64 = (lattitude_deg * cst::DEG_TO_RAD).tan();
    let hour_angle_cos: f64 = -(sun_dec_rad).tan() * tan_lat;
    let hour_angle_rad: f64 = hour_angle_cos.acos();

    let hours_per_radian: f64 = cst::EARTH::SOLAR_DAY / cycle_to_deg * cycle_to_rad;

    let daylight_hours: f64 = hour_angle_rad * hours_per_radian;
    return daylight_hours
}

/// Convert day of year, year to month, day
/// 
/// Inputs
/// ------
/// day_of_year: `u32`
///     Day of year (1-365)
/// 
/// year: `i32`
///     Year (e.g. 2020)
/// 
/// Outputs
/// -------
/// month: `u32`
///     Month (1-12)
/// 
/// day: `u32`
///     Day of month (1-31)
pub fn calc_month_day(
    day_of_year: u32,
    year: i32
) -> (u32, u32) {
    assert!(day_of_year < 366, "Day of year must be less than 366"); 

    let feb_days: u32;
    if check_if_leap_year(year){feb_days = 29;
    } else {feb_days = 28;}

    let month_lengths: Vec<u32> = vec![
        31, feb_days, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    let mut month: u32 = 1;
    let mut sum_days: u32 = month_lengths[0];

    while sum_days < day_of_year -  month_lengths[month as usize - 1]{
        month += 1;
        sum_days += month_lengths[month as usize - 1];
    }

    let month: u32 = month;
    let day: u32 = day_of_year - sum_days;

    return (month, day);
}

/// Convert month, day to day of year
/// 
/// Inputs
/// ------
/// year : `i32`
///     Gregorian Year of common era
/// 
/// month : `i32`
///     Month of year
/// 
/// day : `i32`
///     Day of month
/// 
/// Outputs
/// -------
/// day_of_year: `i32`
///     Day of year
pub fn date_to_doy(
    year: i32, 
    month: i32,
    day: i32
) -> i32 {
    // TODO-TD: add checks to ensure valid date entered

    let mut total_days: i32 = 0;
    let long_months: Vec<i32> = [1, 3, 5, 7, 8, 10, 12].to_vec();
    let short_months: Vec<i32> = [4, 6, 9, 11].to_vec();

    let is_leap_year: bool = check_if_leap_year(year); 

    for i_month in 1..month {
        if long_months.contains(&i_month) {
            total_days += 31;
        } else if short_months.contains(&i_month){
            total_days += 30;
        } else if i_month == 2 && is_leap_year {
            total_days += 29;
        } else {
            total_days += 28;
        }

    }
    return total_days + day
}

/// Check if the year is a leap year
/// 
/// Inputs
/// ------
/// year: `i32`
///     Gregorian Year of common era.
/// 
/// Outputs
/// -------
/// is_leap_year: `bool`
///     True if leap year
fn check_if_leap_year(year: i32) -> bool {
    let rule1: bool = year % 4 == 0;
    let rule2: bool = year % 100 != 0;
    let rule3: bool = year % 400 == 0;
    return rule1 && (rule2 || rule3);
}

/// Convert gregorian date to julian day
/// 
/// Inputs
/// ------
/// year: `i32`
///     Common Era year
/// 
/// month: `i32`
///     Month number of year
/// 
/// day: `i32`
///     Day of month
/// 
/// Outputs
/// -------
/// julian_day_num: `i32`
///     Julian Day Number
pub fn date_to_julian_day_num(
    year: i32,
    month: i32,
    day: i32
) -> i32 {
    // FIXME-TD: fix magic numbers

    let del_month: i32 = (month - 14) / 12; // Adjusts for jul & aug
    let julian_day_num: i32 = (1461 * (year + 4800 + del_month))/4 
        + (367 * (month - 2 - 12 * (del_month)))/12 
        - (3 * ((year + 4900 + del_month) / 100))/4 
        + day - 32075;

    return julian_day_num
}

/// Julian day number to gregorian year, month, day
/// 
/// Inputs
/// ------
/// jd : `i32`
///     Julian day
/// 
/// Outputs
/// -------
/// year : `i32`
///     Gregorian year of common era.
/// 
/// month : `i32`
///     Month of year.
/// 
/// day : `i32`
///     Day of month.
pub fn julian_to_gregorian(
    jd: i32
) -> (i32, i32, i32){
    // FIXME-TD: fix magic numbers
    let e: i32 = (jd + 1401 + (((4*jd + 274277)/146097)* 3) / 4 - 38) * 4 + 3;
    let h: i32 = 5 * (e % 1461) / 4 + 2;

    let day: i32 = (h % 153) / 5 + 1;
    let month: i32 = ((h / 153 + 2) % 12) + 1;
    let year: i32 = e / 1461 - 4716 + (14 - month) / 12;

    return (year, month, day)

}

/// Year Month Day Hour Minute Second to DateTime
/// 
/// Inputs
/// ------
/// year: `i32`
///     Gregorian year of common era.
/// 
/// month: `u32`
///     Month of year.
/// 
/// day: `u32`
///     Day of month.
/// 
/// hour: `u32`
///     Hour of day.
/// 
/// min: `u32`
///     Minute of hour.
/// 
/// sec: `u32`
///     Second of minute.
/// 
/// Outputs
/// -------
/// date_time: `DateTime<Utc>`
///     DateTime object in UTC.
pub fn ymd_hms_to_datetime(
    year: i32,
    month: u32,
    day: u32,
    hour: u32,
    min: u32,
    sec: u32
)-> DateTime<Utc> {
    let date: NaiveDate = NaiveDate::from_ymd_opt(year, month, day).unwrap();
    let time: NaiveTime = NaiveTime::from_hms_opt(hour, min, sec).unwrap();

    let dt: NaiveDateTime = NaiveDateTime::new(date, time);
    let date_time: DateTime::<Utc> = DateTime::from_naive_utc_and_offset(dt, Utc); 
    return date_time
}

/// Calculate schwarzchild radius of a given mass
/// 
/// Inputs
/// ------
/// mass: `f64`
///     Mass of object in kg
/// 
/// Outputs
/// -------
/// radius: `f64`
///     Schwarzchild radius in m
pub fn calc_schwarzchild_radius(
    mass: f64
) -> f64 {
    let radius: f64 = 2.0 * mass * cst::GRAV_CONST / cst::SPEED_OF_LIGHT.powi(2);
    return radius
}

/// Calculate time dilation of relative velocity
/// 
/// Inputs
/// ------
/// t_1: `f64`
///     Time in seconds
/// 
/// rel_vel: `f64`
///     Relative velocity in km/s
/// 
/// Outputs
/// -------
/// t_2: `f64`
///     Dilated time in seconds
pub fn calc_time_dilation(
    t_1: f64,
    rel_vel: f64
) -> f64 {
    let lorentz: f64 = (1.0 - rel_vel.powi(2)/cst::SPEED_OF_LIGHT.powi(2)).sqrt();
    let t_2: f64 = t_1 / (lorentz);
    return t_2
}

/// Calculate apparant angular size of object in field of view
/// 
/// Inputs
/// ------
/// object_radius: `f64`
///     Radius of object in km
/// 
/// radial_distance: `f64`
///     Distance from observer to object in km
pub fn calc_angular_size(
    object_radius: f64,
    radial_distance: f64
) -> f64 {
    let arc_rads: f64 = (object_radius / radial_distance).atan(); 
    return arc_rads
}

/// Rectangular coordinates to geodetic
/// 
/// Inputs
/// ------
/// ecef: `Vector3<f64>
///     Rectangular coordinates in km
/// 
/// Outputs
/// -------
/// lla: `Vector3<f64>`
///     Geodetic coordinates in degrees
pub fn ecef_to_lla(
    ecef: Vector3<f64>,
    body: Body
) -> Vector3<f64> {
    // Zhu's method
    let a: f64 = body.eq_radius;
    let b: f64 = body.eq_radius * (1. - body.eccentricity.powi(2)).sqrt();
    
    let ecc_2: f64 = (a.powi(2) - b.powi(2)) / a.powi(2);
    let ecc_2_prime: f64 = a.powi(2) / b.powi(2) - 1.0;
    
    let x: f64 = ecef[0] * cst::KILO;
    let y: f64 = ecef[1] * cst::KILO;
    let z: f64 = ecef[2] * cst::KILO;

    let p: f64 = (x.powi(2) + y.powi(2)).sqrt();
    let g: f64 = p.powi(2) + (1.0 - ecc_2) * z.powi(2) - 
        ecc_2 * (a.powi(2) - b.powi(2));
    let f: f64 = 54.0 * b.powi(2) * z.powi(2);
    let c: f64 = ecc_2.powi(2) * f * p.powi(2) / (g.powi(3));
    
    let s: f64 = (1.0 + c + (c.powi(2) + 2.0 * c).sqrt()).powf(1. / 3.);
    let cap_p: f64 = f / (3.0 * (s + 1.0 + 1.0 / s).powi(2) * g.powi(2));

    let q: f64 = (1.0 + 2.0 * ecc_2.powi(2) * cap_p).sqrt();
    let r_0: f64 = -cap_p * ecc_2 * p /(1.0+q) + 
        ((a.powi(2)/2.0)*(1.0 + 1.0 / q) - 
        cap_p * (1.0 - ecc_2) * z.powi(2) / (q * (1.0 + q)) - 
        cap_p*p.powi(2)/2.0).sqrt();


    let u: f64 = ((p - ecc_2*r_0).powi(2) + z.powi(2)).sqrt();
    let v: f64 = ((p - ecc_2*r_0).powi(2) + (1.0 - ecc_2)*z.powi(2)).sqrt();
    let z_0: f64 = b.powi(2) * z / (a * v);

    let lat: f64 = cst::to_degrees((z + ecc_2_prime * z_0).atan2(p));
    let lon: f64 =  cst::to_degrees(y.atan2(x));
    let alt: f64 = u * (1.0 - b.powi(2) / (a * v));

    let lla: Vector3<f64> = Vector3::new(lat, lon, alt);
    return lla
}

/// Geodetic to rectangular coordinates
/// E.g. Latitude, Longitude, Altitude to ECEF
/// 
/// Inputs
/// ------
/// lla: `Vector3<f64>`
///     geodetic coords
/// 
/// Outputs
/// -------
/// xyz: `Vector3<f64>`
///     Cartesian coords
pub fn planetodetic_to_cartesian_rotational(
    lla: Vector3<f64>,
    equatorial_radius: f64,
    eccentricity: f64
) -> Vector3<f64> {
    let radius: f64 = calc_prime_vertical(
        lla[0], 
        equatorial_radius,
        eccentricity);
    let x: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].cos();
    let y: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].sin();
    let z: f64 = ((1.0 - eccentricity.powi(2)) * radius + lla[2]) * lla[0].sin();
    let xyz: Vector3<f64> = Vector3::new(x, y, z); 
    return xyz
}

/// Calculate prime vertical radius to surface at latitude
/// 
/// Inputs
/// ------
/// lat_deg: `f64`
///     Lattitude in degrees
pub fn calc_prime_vertical(
    lat_deg: f64, 
    equatorial_radius: f64,
    eccentricity: f64
) -> f64 {
    let lat_radians: f64 = PI * lat_deg / 180.0;
    let radius: f64 = 
        equatorial_radius / (1.0 - (eccentricity * lat_radians.sin()).powi(2)).sqrt();
    return radius
}

/// Map between fixed frame observation to enu
/// 
/// Inputs
/// ------
/// pos_lla: `Vector3<f64>`
///     Lattitude, Longitude, Altitude
/// 
/// ecef_2: `Vector3<f64>`
///     ECEF object
/// 
/// Outputs
/// -------
/// enu: `Vector3<f64>`
///     East, North, Up
pub fn fixed_frame_to_enu(
    observer_lla: Vector3<f64>, 
    p_tgt_fixed: Vector3<f64>,
    equatorial_radius: f64,
    eccentricity: f64
) -> Vector3<f64> {
    let observer_ecef: Vector3<f64> = planetodetic_to_cartesian_rotational(
        observer_lla,
        equatorial_radius,
        eccentricity);
    let vec_ecef: Vector3<f64> = p_tgt_fixed - observer_ecef;
    let ecef_enu: Matrix3<f64> = Matrix3::new(
        -observer_lla[1].sin(), observer_lla[1].cos(), 0.0,
        -observer_lla[1].cos()*observer_lla[0].sin(), -observer_lla[1].sin()*observer_lla[0].sin(), observer_lla[0].cos(),
        observer_lla[1].cos()*observer_lla[0].cos(), observer_lla[1].sin()*observer_lla[0].cos(), observer_lla[0].sin());
    let enu: Vector3<f64> = ecef_enu * vec_ecef;
    return enu
}

/// Map between enu and fixed frame
/// 
/// Inputs
/// ------
/// pos_lla: `Vector3<f64>`
///     Lattitude, Longitude, Altitude
/// enu: `Vector3<f64>`
///     East, North, Up
/// 
/// Outputs
/// -------
/// ecef: `Vector3<f64>`
///     ECE
pub fn enu_to_ecef(
    pos_lla: Vector3<f64>, 
    enu: Vector3<f64>,
    equatorial_radius: f64,
    eccentricity: f64
) -> Vector3<f64> {
    let enu_ecef: Matrix3<f64> = Matrix3::new(
        -pos_lla[1].sin(), -pos_lla[1].cos()*pos_lla[0].sin(), pos_lla[1].cos()*pos_lla[0].cos(),
        pos_lla[1].cos(), -pos_lla[1].sin()*pos_lla[0].sin(), pos_lla[1].sin()*pos_lla[0].cos(),
        0.0, pos_lla[0].cos(), pos_lla[0].sin()
    );
    let vec_ecef: Vector3<f64> = enu_ecef * enu;
    let pos_ecef: Vector3<f64> = planetodetic_to_cartesian_rotational(
        pos_lla,
        equatorial_radius,
        eccentricity);
    let ecef: Vector3<f64> = vec_ecef - pos_ecef;
    return ecef
}

#[cfg(test)]
mod spacetime_tests {
    use super::*;

    #[test]
    fn test_date_to_doy(){
        let year: i32 = 2023;
        let month: i32 = 2;
        let day: i32 = 5;

        let doy: i32 = date_to_doy(year, month, day);
        assert_eq!(doy, 36);
    }

    #[test]
    fn test_julian_day(){
        let year: i32 = 2023;
        let month: i32 = 2;
        let day: i32 = 5;

        let julian_day: i32 = date_to_julian_day_num(year, month, day);
        assert_eq!(julian_day, 2459981);
    
    }

    #[test]
    fn test_earth_day_length(){
        let year: i32 = 2023;
        let month: i32 = 3;
        let day: i32 = 5;
        let long_deg: f64 = -83.7101;
        let lat_deg: f64 = 42.2929; 

        let julian_day: i32 = date_to_julian_day_num(year, month, day);
        let day_light_hrs: f64 = calc_earth_day_length(lat_deg, long_deg, julian_day);
        assert_eq!(day_light_hrs, 8.403851475586468);

    }


    #[test]
    fn test_schwarzchild_radius(){
        let mass: f64 = 1.0;
        let radius: f64 = calc_schwarzchild_radius(mass);
        assert_eq!(radius, 1.48523238761875e-27);
    }


    #[test]
    fn test_ecef_to_lla(){
        let pos_ecef = Vector3::new(
            -2383.0,
            -4662.0,
            5124.0
        );
        let earth: Body = Body {
            name: String::from("Earth"),
            grav_param: cst::EARTH::GRAV_PARAM,
            eq_radius: cst::EARTH::RADIUS_EQUATOR,
            rotation_rate: cst::EARTH::ROT_RATE,
            sidereal_day_hours: cst::EARTH::SIDEREAL_DAY,
            eccentricity: cst::EARTH::SURFACE_ECC
        };

        let pos_lla = ecef_to_lla(
            pos_ecef, earth);

        assert_eq!(pos_lla, Vector3::new(
            44.549281984923006, 
            -117.07402908139512, 
            958212.9075099035));
    }

    #[test]
    fn test_lla_to_ecef(){
        let pos_lla = Vector3::new(
            44.54968779193849, 
            -117.07402908139512, 
            958238.4011472173);
        let pos_ecef = planetodetic_to_cartesian_rotational(
            pos_lla, 
            cst::EARTH::RADIUS_EQUATOR, 
            cst::EARTH::SURFACE_ECC);
        assert_eq!(pos_ecef, Vector3::new(
            -4157947.6438792264, 
            4593265.390986962, 
            3925488.3377592554));
    }
}
