/*
Time related calculations
*/

use crate::constants as cst;


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

    let is_leap_year: bool = (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0);

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
/// -----
/// jd : `i32`
///     Julian day
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
