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



/// Stumpff function
/// 
/// Parameters
/// ----------
/// x:
/// 
/// order:
/// 
/// eps: 
pub fn stumpff(x: f64, order: i32, eps:  f64) ->  f64 {
    (-1.).powi(n)*x.powi(n) / (order + 2*n).factorial();
};

pub fn stumpff_k(
    t: f64,
    k: f64,
    N: i32=1,
    delta: f64=1e-15
) -> f64 {
    let mut sk:  f64;
    while abs(t**N/factorial(2*N+k+2))>delta:N+=5{
        if n<N{
            sk=lambda n:t/((2*n+k+1)*(2*n+k+2))*(1-sk(n+1)); 
        else{
            sk=0.;
        }

    }
    return (1-sk(0))/factorial(k)
}

// def ck_noerror(t,k,N=20):
//     sk=lambda n:t/((2*n+k+1)*(2*n+k+2))*(1-sk(n+1)) if n<N else 0
//     return (1-sk(0))/factorial(k)

// def error_ck(t,k,N):
//     return abs(t**N/factorial(2*N+k+2))

fn ck_analytic(t:  f64) -> Vec<f64> {
    if t==0.{
        return vec![1., 1., 0.5, 1./6.]
    } else{

        let y = t.abs().powf(0.5);
        if t>=0. {
            let c0 = y.cos();
            let c1 = y.sin()/y;
            let c2 = (1 - y.cos())/y**2;
            let c3 = (y - y.sin())/y**3;
            return vec![c0,c1,c2,c3]
        }else {
            let c0 = cosh(y);
            let c1 = sinh(y)/y;
            let c2 = -(1-cosh(y))/y**2;
            let c3 = -(y-sinh(y))/y**3;
            return vec![c0,c1,c2,c3]
        }
        
        }
    }
