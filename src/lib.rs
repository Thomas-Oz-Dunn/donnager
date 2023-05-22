
pub mod donnager;


use pyo3::prelude::*;

#[pyfunction]
fn calc_earth_day_length(: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pymodule]
fn donnager(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calc_earth_day_length, m)?)?;
    Ok(())
}