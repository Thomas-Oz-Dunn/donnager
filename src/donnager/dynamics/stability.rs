// Stability Analysis

use nalgebra::*;

/// Calculate lyanpunov stablity
/// 
/// Inputs
/// ------
/// 
pub fn calc_lyanpunov_stability(
    mass_matrix: MatrixN<f64>,
    damp_matrix: MatrixN<f64>,
    gyro_matrix: MatrixN<f64>,
    stiff_matrix: MatrixN<f64>
) -> Vec<Tuple<f64, f64>> {
    let mass_inv = mass_matrix.try_inverse().unwrap();

    let a_mat = MatrixN::new(
        [MatrixN::identity(), MatrixN::zeros()],
        [[-mass_inv * stiff_matrix],[-mass_matrix * (damp_matrix + gyro_matrix)]]
    );
    // Eigenvalue decomposition
    
}