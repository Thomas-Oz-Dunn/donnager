

#[derive(Clone, Debug, PartialEq)]
pub struct Vehicle{
    pub name: String,
    pub mass_0: f64,
    pub mass_prop: f64,
    pub engine_type: String,
    pub engine_isp: f64,
}

pub struct Stage{
    pub vehicle: Vehicle,
    pub struct_coeff: f64
}


pub struct Multistage{
    pub name: String,
    pub n_stage: f64,
    pub glow: f64,
    pub stages: [Stage]
}




impl Multistage {

}