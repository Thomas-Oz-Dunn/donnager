/*
Space time related functions
 */

use nalgebra as na;
use na::{Vector3, Matrix3};
use std::f64::consts::PI;
use chrono::{DateTime, Utc, Datelike, Timelike};

use crate::constants as cst;

/// Gravitational Body
#[derive(Clone, Debug, PartialEq)]
pub struct Body{
    pub name: String,
    pub grav_param: f64,
    pub eq_radius: f64,
    pub rotation_rate: f64,
    pub eccentricity: f64
}

pub enum ReferenceFrames {
    ECI,
    ECEF,
    PFCL,
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
        let particle: Particle = Particle {mass: mass, motion: motion};
        particle
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
        let vel: f64 = (2.0 * self.grav_param / radius).sqrt();
        return vel
    }

    /// Calculate period of orbit at semi major axis
    /// 
    /// Inputs
    /// ------
    /// semi_major_axis: `f64`
    ///     Semi major axis of orbital ellipse
    /// 
    /// Outputs
    /// -------
    /// period: `f64`
    ///    Period of orbit in seconds
    pub fn calc_period(&self, semi_major_axis: f64) -> f64 {
        let period: f64 = 2.0 * PI * (semi_major_axis.powi(3)/self.grav_param).sqrt();
        return period
    }
    
    /// Calculate radius for stationary orbit above body surface
    /// 
    /// Outputs
    /// -------
    /// r_mag: `f64`
    ///     Magnitude of radius for stationary orbit
    pub fn calc_stationary_orbit(&self) -> f64 {
        let period: f64 = 2.0 * PI / self.rotation_rate;
        let a: f64 = self.grav_param * period.powi(2); 
        let r_mag: f64 = (a / (4.0 * PI.powi(2))).powf(1.0 / 3.0); 
        return r_mag
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
    pub fn geodetic_to_xyz(&self, lla: Vector3<f64>) -> Vector3<f64> {
        let radius: f64 = self.calc_prime_vertical(lla[0]);
        let x: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].cos();
        let y: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].sin();
        let z: f64 = ((1.0 - self.eccentricity.powi(2)) * radius + lla[2]) * lla[0].sin();
        let xyz: Vector3<f64> = Vector3::new(x, y, z); 
        return xyz
    }

    /// Calculate prime vertical radius to surface at latitude
    /// 
    /// Inputs
    /// ------
    /// lat_deg: `f64`
    ///     Lattitude in degrees
    pub fn calc_prime_vertical(&self, lat_deg: f64) -> f64 {
        let lat_radians: f64 = PI * lat_deg / 180.0;
        let radius: f64 = 
            self.eq_radius / (1.0 - (self.eccentricity * lat_radians.sin()).powi(2)).sqrt();
        return radius
    }

    /// Rectangular coordinates to geodetic
    /// E.g. ECEF to LLH
    pub fn xyz_to_geodetic(&self, xyz: Vector3<f64>) -> Vector3<f64> {
        // Zhu's method
        let a: f64 = self.eq_radius;
        let ecc_2: f64 = self.eccentricity.powi(2);

        let b: f64 = (a.powi(2)*(1.0 - ecc_2)).sqrt();
        let ecc_2_prime: f64 = a.powi(2) / b.powi(2) - 1.0;
        let p: f64 = (xyz[0].powi(2) + xyz[1].powi(2)).sqrt();
        let g: f64 = p.powi(2) + (1.0 - ecc_2) * xyz[2].powi(2) - 
            ecc_2 * (a.powi(2) - b.powi(2));
        let c: f64 = ecc_2.powi(2) * 54.0 * b.powi(2) * xyz[2].powi(2) * p.powi(2) / (g.powi(3));
        let s: f64 = (1.0 + c + (c.powi(2) + 2.0 * c).sqrt()).powf(1.0 / 3.0);
        let cap_p: f64 = 54.0 * b.powi(2)*xyz[2].powi(2) / 
            (3.0 * (s + 1.0 + 1.0 / s).powi(2) * g.powi(2));
        let q: f64 = (1.0 + 2.0 * ecc_2.powi(2) * cap_p).sqrt();
        let r_0: f64 = -cap_p * ecc_2 * p /(1.0+q) + 
            ((a.powi(2)/2.0)*(1.0 + 1.0 / q) - 
            cap_p * (1.0 - ecc_2) * xyz[2].powi(2) / (q * (1.0 + q)) - 
            cap_p*p.powi(2)/2.0).sqrt();
        let u: f64 = ((p - ecc_2*r_0).powi(2) + xyz[2].powi(2)).sqrt();
        let v: f64 = ((p - ecc_2*r_0).powi(2) + (1.0 - ecc_2)*xyz[2].powi(2)).sqrt();
        let z_0: f64 = b.powi(2) * xyz[2] / (a * v);

        let alt: f64 = u * (1.0 - b.powi(2) / (a * v));
        let lat: f64 = ((xyz[2] + ecc_2_prime*z_0)/p).atan();
        let lon: f64 = (xyz[1] / xyz[0]).atan();
        let lla: Vector3<f64> = Vector3::new(lat, lon, alt);
        return lla
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
    pub fn to_body(&self, name: String, eq_radius: f64, rotation_rate: f64, eccentricity: f64) -> Body {
        let grav_param: f64 = self.mass * cst::GRAV_CONST;
        let body: Body = Body {
            name,
            grav_param,
            eq_radius,
            rotation_rate,
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


#[derive(Clone, Debug, PartialEq)]
pub enum Frame{
    ENU,            // ENU
    ECEF,        // ECEF
    ECI,     // ECI
    LlaE,               // LLA (Earth)
    MCMF,       // MCMF
    MCI,    // MCI
    LlaM,              // LLA (Mars)
    Heliocentric            // Sun
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
        let prime_vertical: f64 = self.body.calc_prime_vertical(self.pos_lla[0]);
        let dir: Vector3<f64> = self.body.geodetic_to_xyz(self.pos_lla);
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


    /// Map between fixed frame observation to enu
    pub fn ecef_to_enu(&self, ecef_2: Vector3<f64>) -> Vector3<f64> {
        let pos_lla: Vector3<f64> = self.pos_lla;
        let pos_ecef: Vector3<f64> = self.body.geodetic_to_xyz(pos_lla);
        let vec_ecef: Vector3<f64> = ecef_2 - pos_ecef;
        let ecef_enu: Matrix3<f64> = Matrix3::new(
            -pos_lla[1].sin(), pos_lla[1].cos(), 0.0,
            -pos_lla[1].cos()*pos_lla[0].sin(), -pos_lla[1].sin()*pos_lla[0].sin(), pos_lla[0].cos(),
            pos_lla[1].cos()*pos_lla[0].cos(), pos_lla[1].sin()*pos_lla[0].cos(), pos_lla[0].sin());
        let enu: Vector3<f64> = ecef_enu * vec_ecef;
        return enu
    }

    /// Map between enu and fixed frame
    /// 
    /// 
    pub fn enu_to_ecef(&self, enu: Vector3<f64>) -> Vector3<f64> {
        let pos_lla: Vector3<f64> = self.pos_lla;
        let enu_ecef: Matrix3<f64> = Matrix3::new(
            -pos_lla[1].sin(), -pos_lla[1].cos()*pos_lla[0].sin(), pos_lla[1].cos()*pos_lla[0].cos(),
            pos_lla[1].cos(), -pos_lla[1].sin()*pos_lla[0].sin(), pos_lla[1].sin()*pos_lla[0].cos(),
            0.0, pos_lla[0].cos(), pos_lla[0].sin()
        );
        let vec_ecef: Vector3<f64> = enu_ecef * enu;
        let pos_ecef: Vector3<f64> = self.body.geodetic_to_xyz(self.pos_lla);
        let ecef: Vector3<f64> = vec_ecef - pos_ecef;
        return ecef
    }

}

// Calculate eci to ecef rotation matrix
/// 
/// Inputs
/// ------
/// date_time
pub fn calc_eci_ecef_rotam(date_time: DateTime<Utc>) -> Matrix3<f64> {
    let year = date_time.year();
    let month = date_time.month() as i32;
    let day = date_time.day() as i32;

    let hours = date_time.hour() as f64;
    let minutes = date_time.minute() as f64;
    let seconds = date_time.second() as f64;

    let julian_day: f64 = date_to_julian_day_num(year, month, day) as f64;

    let sidereal_time: f64 = (hours + (minutes + seconds / 60.) / 60.) / 24.;

    let rot_rate_rad_day: f64 = 
        cst::EARTH_ROT_RATE * 3600. * cst::EARTH_SIDEREAL_DAY; 
    let theta: f64 = 
        rot_rate_rad_day * (julian_day + sidereal_time - cst::J2000_DAY);

    let rotam: Matrix3<f64> = Matrix3::<f64>::new(
        theta.cos(), -theta.sin(), 0.,
        theta.sin(), theta.cos(), 0.,
        0., 0., 1.);

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
pub fn calc_day_length(
    lat_deg: f64,
    long_deg: f64,
    julian_day: i32
) -> f64 {
    let j2000_days: i32 = julian_day - 2451545;
    
    let j_star: f64 = (j2000_days as f64) - long_deg / 360.;
    let m: f64 = (357.5291 + 0.98560028 * j_star) % 360.;
    let eq_cen: f64 = 1.9148*(m.sin()) + 0.02*((2.*m).sin()) + 0.0003*((3.*m).sin());
    let lambda: f64 = (m + eq_cen + 180. + cst::EARTH_ARG_PERIHELION) % 360.;
    let sun_dec_rad: f64 = (lambda.sin() * (cst::EARTH_AXIAL_TILT).sin()).asin();
    let tan_lat: f64 = (lat_deg * cst::DEG_TO_RAD).tan();

    let hour_angle_rad: f64 = (-(sun_dec_rad).tan() * tan_lat).acos();
    let hours: f64 = hour_angle_rad / 15.0 * cst::RAD_TO_DEG;
    return hours
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

    while sum_days < day_of_year {
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
pub fn date_to_julian_day_num(
    year: i32,
    month: i32,
    day: i32
) -> i32 {

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
    let e: i32 = (jd + 1401 + (((4*jd + 274277)/146097)* 3) / 4 - 38) * 4 + 3;
    let h: i32 = 5 * (e % 1461) / 4 + 2;

    let day: i32 = (h % 153) / 5 + 1;
    let month: i32 = ((h / 153 + 2) % 12) + 1;
    let year: i32 = e / 1461 - 4716 + (14 - month) / 12;

    return (year, month, day)

}


#[cfg(test)]
mod time_tests {
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
    fn test_day_length(){
        let year: i32 = 2023;
        let month: i32 = 2;
        let day: i32 = 5;
        let long_deg: f64 = -83.7101;
        let lat_deg: f64 = 42.2929; 

        let julian_day: i32 = date_to_julian_day_num(year, month, day);
        let day_light_hrs: f64 = calc_day_length(lat_deg, long_deg, julian_day);
        assert_eq!(day_light_hrs, 8.54933135165009);

    }
}
