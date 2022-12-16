// 

pub fn calc_day_length(
    lattitude: f64,
    doy: i32
) -> f64 {
    let a: f64 = (0.98565 * ((doy - 2) as f64)).sin();
    let b: f64 = 360 / PI * 0.0167 * a;
    let c: f64 = 360 / 365.24 * ((doy + 10) as f64);
    let sun_declination: ((b + c).cos() * (-23.44).sin()).asin();
    let hour_angle: f64 = (-sun_declination.tan() * lattitude.tan()).acos();
    let hours = hour_angle / 15.0;
    return hours
}