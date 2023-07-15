/// Hypergeometric function 2F1(3, 1, 5/2, x), see [Battin].
/// https://en.wikipedia.org/wiki/Hypergeometric_function
pub fn hyp2f1b(x: f64) -> f64 {

    if x >= 1.0{
        return f64::INFINITY;
    }else{
        let mut res: f64 = 1.0;
        let mut res_old: f64;
        let mut term: f64 = 1.0;
        let mut ii: f64 = 0.;

        while term > 0.{
            term = term * (3. + ii) * (1. + ii) / (5. / 2. + ii) * x / (ii + 1.);
            res_old = res;
            res += term;
            ii += 1.;
        }
        return res

    }
}

/// Analytic Stumpff function 
/// 
/// Parameters
/// ----------
/// t: f64
///     
/// Returns
/// -------
/// c_n: Vec<f64>
///     Four orders of solutions
pub fn stumpff_analytic(t:  f64) -> Vec<f64> {
    if t==0.{
        return vec![1., 1., 0.5, 1./6.]
    } else{
        let y = t.abs().powf(0.5);
        if t>=0. {
            let c_0: f64 = y.cos();
            let c_1: f64 = y.sin()/y;
            let c_2: f64 = (1. - y.cos())/y.powi(2);
            let c_3: f64 = (y - y.sin())/y.powi(3);
            return vec![c_0, c_1, c_2, c_3]
        }else{
            let c_0: f64 = y.cosh();
            let c_1: f64 = y.sinh()/y;
            let c_2: f64 = -(1. - y.cosh())/y.powi(2);
            let c_3: f64 = -(y - y.sinh())/y.powi(3);
            return vec![c_0, c_1, c_2, c_3]
        }
    }
}


fn error_ck(t: f64, k:  u32,  N:  u32) -> f64{
    let base: u32 = 2 * N + k  + 2;
    let eps = (t.powi(N as i32) / (factorial(base) as f64)).abs() ;
    return eps
}

/// Factorial
/// 
/// Calculate n!
/// 
/// Inputs
/// ------
/// n: u32
///     Base to factorialize 
pub fn factorial(n: u32) -> u32 {
    return (1..=n).product()
}

#[cfg(test)]
mod math_tests {

    use super::*;

    #[test]
    fn test_stumpff_analytic(){
        let t_n1 = -1.0;
        let c_n_n1 = stumpff_analytic(t_n1);

        let t_0 = 0.0;
        let c_n_0 = stumpff_analytic(t_0);

        let t_1= 1.0;
        let c_n_1 = stumpff_analytic(t_1);
    }

    #[test]
    fn test_factorial(){
        let n_1: u32 = 1;
        let n_1_fact = factorial(n_1);

        assert_eq!(1, n_1_fact);

        let n_10: u32 = 10;
        let n_fact_10 = factorial(n_10);

        assert_eq!(3628800, n_fact_10);
    }
}