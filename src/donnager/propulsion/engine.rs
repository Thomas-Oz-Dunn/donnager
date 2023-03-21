
#[derive(Clone, Debug, PartialEq)]
pub struct Engine{
    pub name: String,
    pub engine_type: EngineType,
    pub isp: f64
}

#[derive(Clone, Debug, PartialEq)]
pub enum EngineType{
    Electric,
    Chemical,
    Nuclear
}

impl Engine {
    
}