use argmin::core::{CostFunction, Error};
use argmin::solver::brent::BrentRoot;
use argmin::core::Executor;

pub struct GeneralCost<F>
where
    F: Fn(f64) -> f64,
{
    pub func: F,
}

impl<F> CostFunction for GeneralCost<F>
where
    F: Fn(f64) -> f64,
{
    type Param = f64;
    type Output = f64;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        let result = (self.func)(*p);
        Ok(result)
    }
}

pub fn find_interval<F>(dens_guess: f64, f: F) -> (f64, f64)
where 
    F: Fn(f64) -> f64,
{
    let factor = 1.5;
    let ntry = 50;

    let mut a = dens_guess*0.8;
    let mut b = dens_guess*1.2;
    for _ in 0..ntry{
        let fa = f(a);
        let fb = f(b);

        if fa * fb <= 0. {
            return (a, b);
        }
        if fa.abs() < fb.abs() {
            a += factor * (a - b);
        } else {
            b += factor * (b - a);
        }
    }
    (0., 0.)
}

pub fn find_root<F>(interval: (f64, f64), f: F) -> f64 
where 
    F: Fn(f64) -> f64,
{
    let cost = GeneralCost { func: f };
    let solver = BrentRoot::new(interval.0, interval.1, 1e-11);
    let res = Executor::new(cost, solver)
        .configure(|state| state.param(interval.0).max_iters(100))
        .run()
        .unwrap();
    let solution = match res.state.param {
        Some(p) => p,
        None => 0.,
    };
    solution
}