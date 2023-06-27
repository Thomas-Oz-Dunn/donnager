
pub mod donnager;

// use pyo3::prelude::*;

// #[pyfunction]
// fn calc_earth_day_length(
//     lattitude_deg: f64, 
//     longitude_deg: f64, 
//     julian_day: i32
// ) -> PyResult<f64> {
//     let hours = donnager::spacetime::calc_earth_day_length(
//         lattitude_deg, 
//         longitude_deg, 
//         julian_day);
//     return PyResult(hours)
// }

// #[pymodule]
// fn donnager(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(calc_earth_day_length, m)?)?;
//     Ok(())
// }