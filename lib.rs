

use pyo3::prelude::*;

use donnager as donnager_rs;

#[pyfunction]
fn calc_earth_day_length(
    lattitude_deg: f64, 
    longitude_deg: f64, 
    julian_day: i32
) -> PyResult<f64> {
    let hours = donnager_rs::spacetime::calc_earth_day_length(
        lattitude_deg, 
        longitude_deg, 
        julian_day);
    return Ok(hours)
}

#[pymodule]
fn donnagerpy(
    _py: Python<'_>, 
    m: &PyModule
)-> PyResult<()> {
    let spacetime = PyModule::new(_py, "submodule")?;
    spacetime.add_function(wrap_pyfunction!(calc_earth_day_length, m)?)?;
    m.add_submodule(spacetime)?;
    Ok(())
}
