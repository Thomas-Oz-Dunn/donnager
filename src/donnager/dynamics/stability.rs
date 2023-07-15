// Stability Analysis

use nalgebra::*;

/// Calculate lyanpunov stablity
/// 
/// Inputs
/// ------
/// mass_matrix
/// 
/// damp_matrix
/// 
/// gyro_matrix
/// 
/// stiff_matrix
pub fn calc_2d_lyanpunov_stability(
    mass_matrix: Matrix2<f64>,
    damp_matrix: Matrix2<f64>,
    gyro_matrix: Matrix2<f64>,
    stiff_matrix: Matrix2<f64>
) -> Vector4<nalgebra::Complex<f64>> {
    let mass_inv: Matrix2<f64> = mass_matrix.try_inverse().unwrap(); 
    let iden: Matrix2<f64> = Matrix2::<f64>::identity(); 
    let zer: Matrix2<f64> = Matrix2::<f64>::zeros(); 
    let left_bot: Matrix2<f64> = -mass_inv * stiff_matrix; 
    let right_bot: Matrix2<f64> = -mass_inv * (damp_matrix + gyro_matrix); 

    let a_mat: Matrix4<f64> = Matrix4::new(
        iden.m11, iden.m12, zer.m11, zer.m12,
        iden.m21, iden.m22, zer.m21, zer.m22,
        left_bot.m11, left_bot.m12, right_bot.m11, right_bot.m12,
        left_bot.m21, left_bot.m22, right_bot.m21, right_bot.m22
    );

    // Complex eigenvalue decomposition
    let eig: Vector4<nalgebra::Complex<f64>> = a_mat.complex_eigenvalues();
    return eig;
}



#[cfg(test)]
mod stability_tests {

    use super::*;
    use crate::donnager::constants as cst;    

    #[test]
    fn test_calc_lyanpunov_stability()
    {
        let mass_ratio = cst::EARTH::MASS / (cst::EARTH::MASS + cst::SUN::MASS);

        // Libration point
        let mass_matrix = Matrix2::<f64>::identity();
        let damp_matrix = Matrix2::<f64>::zeros();
        let gyro_matrix = Matrix2::<f64>::new(0., -2., 2., 0.);
        let stiff_matrix = Matrix2::<f64>::new(
            -3./4., -3.* (3f64.sqrt())/2. * (mass_ratio - 0.5),
            3.* (3f64.sqrt())/2. * (mass_ratio - 0.5), -9./4.
        );

        let roots = calc_2d_lyanpunov_stability(
            mass_matrix,
            damp_matrix,
            gyro_matrix,
            stiff_matrix
        );

    }
}
