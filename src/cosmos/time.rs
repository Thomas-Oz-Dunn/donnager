
use std::f64::consts::PI;

// Calculate hours of sunlight at lattitude
pub fn calc_day_length(
    lattitude: f64,
    doy: i32
) -> f64 {
    let a: f64 = (0.98565 * ((doy - 2) as f64)).sin();
    let c: f64 = 360.0 / 365.24 * ((doy + 10) as f64) + 360.0 / PI * 0.0167 * a;
    let sun_declination: f64 = (c.cos() * (-23.44 as f64).sin()).asin();
    let hour_angle: f64 = (-sun_declination.tan() * lattitude.tan()).acos();
    let hours = hour_angle / 15.0;
    return hours
}

// Convert month, day to day of year
pub fn convert_date_to_doy(
    month: i32,
    day: i32
) -> i32 {
    let mut total_days: i32 = 0;
    let long_months: Vec<i32> = [1, 3, 5, 7, 8, 10, 12].to_vec();
    let short_months: Vec<i32> = [4, 6, 9, 11].to_vec();
    for i_month in 1..month {
        if long_months.contains(&i_month) {
            total_days += 31;
        } else if short_months.contains(&i_month){
            total_days += 30;
        } else {
            // Leap year check?
            total_days += 28;
        }
    }
    return total_days + day
}