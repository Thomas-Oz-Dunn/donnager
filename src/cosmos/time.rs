/*
Time related calculations
*/

use crate::constants as cst;


/// Calculate hours of sunlight at lattitude
/// 
/// Inputs
/// ------
/// lat_deg : `f64`
///     Lat of coordinates in degrees
/// 
/// doy : `i32`
///     Day of year
pub fn calc_day_length(
    lat_deg: f64,
    doy: i32
) -> f64 {
    // Spatial contribution
    let lat_rad: f64 = lat_deg * cst::DEG_TO_RAD;

    // Temporal contribution
    let alpha: f64 = (0.98565 * ((doy - 2) as f64)).sin();
    let alpha_rad: f64 = 0.0334 * alpha * cst::DEG_TO_RAD;
    let solstice_angle: f64 = (doy + 10) as f64 / cst::EARTH_DAYS_PER_YEAR;
    let gamma: f64 = 360.0 * solstice_angle + alpha_rad;
    let sun_dec: f64 = (gamma.cos() * (cst::EARTH_AXIAL_TILT).sin()).asin();
    
    let hour_angle: f64 = ((-sun_dec).tan() * (lat_rad).tan()).acos();
    let hours: f64 = hour_angle / 15.0;
    return hours
}

/// Convert month, day to day of year
/// 
/// Inputs
/// ------
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

/// Convert gregorian datetime to julian time
/// 
/// Inputs
/// ------
/// gregorian_year: `i32`
///     Common Era year
/// 
/// gregorian_month: `i32`
///     Month number of year
/// 
/// day: `i32`
///     Day of month
pub fn gregorian_date_to_julian_day_num(
    gregorian_year: i32,
    gregorian_month: i32,
    day: i32
) -> i32 {

    let del_month: i32 = (gregorian_month - 14) / 12; // Adjusts for jul & aug
    let julian_day_num: i32 = (1461 * (gregorian_year + 4800 + del_month))/4 
        + (367 * (gregorian_month - 2 - 12 * (del_month)))/12 
        - (3 * ((gregorian_year + 4900 + del_month) / 100))/4 
        + day - 32075;

    return julian_day_num
}


// /// Julian day number to gregorian year, month, day
// pub fn julian_to_gregorian(
//     julian_day: i32
// ) -> (i32, i32, i32){

//     return None
// }




#[cfg(test)]
mod time_tests {
    use super::*;

    #[test]
    fn test_date_to_doy(){
        let year = 2023;
        let month = 2;
        let day = 5;

        let doy = date_to_doy(year, month, day);
        assert_eq!(doy, 36);
    }
}