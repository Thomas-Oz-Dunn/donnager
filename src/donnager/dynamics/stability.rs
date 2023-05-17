// Stability Analysis

use nalgebra::*;

use crate::donnager::constants as cst;

/// Calculate lyanpunov stablity
/// 
/// Inputs
/// ------
/// 
pub fn calc_2d_lyanpunov_stability(
    mass_matrix: Matrix2<f64>,
    damp_matrix: Matrix2<f64>,
    gyro_matrix: Matrix2<f64>,
    stiff_matrix: Matrix2<f64>
) -> Vec<f64> {
    let mass_inv = mass_matrix.try_inverse().unwrap(); 
    let iden = Matrix2::<f64>::identity(); 
    let zer = Matrix2::<f64>::zeros(); 
    // Eigenvalue decomposition
    let eig = a_mat.complex_eigenvalues();
    return eig;
}



#[cfg(test)]
mod stability_tests {

    use super::*;

    #[test]
    fn test_calc_lyanpunov_stability()
    {
        let mass_ratio = cst::EARTH::MASS / (cst::EARTH::MASS + cst::SUN::MASS);

        let mass_matrix = Matrix2::<f64>::identity();
        let damp_matrix = Matrix2::<f64>::zeros();
        let gyro_matrix = Matrix2::<f64>::new(0, -2, 2, 0);
        let stiff_matrix = Matrix2::<f64>::new(
            -3./4., -3.* (3f64.powf(0.5))/2. * (mass_ratio - 0.5),
            3.* (3f64.powf(0.5))/2. * (mass_ratio - 0.5), -9./4.
        );

        let roots = calc_lyanpunov_stability(
            mass_matrix,
            damp_matrix,
            gyro_matrix,
            stiff_matrix
        );

    }
}
