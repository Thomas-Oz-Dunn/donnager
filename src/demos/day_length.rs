/*
Hours of sunlight based on lattitude and day
*/

use nalgebra::Vector3;

use donnager::donnager::spacetime as xyzt;

fn main() {
    // Create mesh
    let days = 0..365;
    let lattitudes = -90..90;

    // Calculate hours of sunlight
    let mut sun_hours = Vec::new();

    days.iter().map(|day| {
        lattitudes.iter().map(|lat| {
            sun_hours.push(xyzt::calc_day_length(day, lat))})
        }
    )

    // Plot results

}
