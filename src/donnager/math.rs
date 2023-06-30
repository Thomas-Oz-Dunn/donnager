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

        while(term > 0.){
            term = term * (3. + ii) * (1. + ii) / (5. / 2. + ii) * x / (ii + 1.);
            res_old = res;
            res += term;
            ii += 1.;
        }
        return res

    }
}
