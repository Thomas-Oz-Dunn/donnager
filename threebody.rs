/*
n-body simulation demonstration
Barnes Hut
*/

pub struct Particle {
    pub pos: Vec3<f64>,
    pub vel: Vec3<f64>,
    pub mass: f64,
}

pub struct Node {
    pub mass: f64
}

impl Node {

}

pub struct Tree {
    pub nodes: Vec<Node>,
    pub center_of_mass: Vec3<f64>
}

impl Tree {

}


/// Calculate 2-body acceleration, brute force
fn calc_acceleration(
    particle1: Particle,
    particle2: Particle,
    field_strength: f64
) -> (Vec3<f64>, Vec3<f64>) {
    let distance: Vec3<f64> = particle1.pos - particle2.pos;
    let radius: f64 = distance.norm();
    let force: f64 =  field_strength * (particle1.mass + particle2.mass) / radius.powi(2);
    let dir: Vec3<f64> = distance / radius;
    let acc_1: Vec3<f64> = -force / particle1.mass * dir;
    let acc_2: Vec3<f64> = force / particle2.mass * dir;
    return (acc_1, acc_2)
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
    let t_end: u_size = 60;
    let field_strength: f64 = 4.0;

    // Control duration and precision
    for dt in (0..t_end).iter(){
        (acc1, acc2) = calc_acceleration(particle1, particle2, field_strength);

        particle1.pos += particle1.vel * dt;
        particle2.pos += particle2.vel * dt;

        particle1.vel += acc1 * dt;
        particle2.vel += acc2 * dt;

    }

    // 3 body

}