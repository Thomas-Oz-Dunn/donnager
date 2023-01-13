/*
n-body simulation demonstration
Barnes Hut
*/

pub struct Particle {
    pub pos: Vec3<f64>,
    pub vel: Vec3<f64>,
    pub mass: f64,
}


// NOTE: Make agnostic to 1, 2, 3, N d space?

pub struct Node {}

impl Node {
    
}

pub struct Tree {
    pub nodes: Vec<Node>
}

impl Tree {

}







fn main() {

    // 1 body
    let mut particle1: Particle = Particle {
        pos = Vec3::new(0.,0.,0.),
        vel = Vec3::new(0.,0.,1.),
        mass = 5.0
    };

    // 2 body
    let mut particle2: Particle = Particle {
        pos = Vec3::new(1.,0.,0.),
        vel = Vec3::new(0.,-2.,0.),
        mass = 5.0
    };

    
    // Run


    // 3 body

}