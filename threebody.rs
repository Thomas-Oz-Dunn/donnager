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


/// Calculate acceleration, brute force
fn calc_acceleration(
    particles: Vec<Particle>,
    field_strength: f64
) -> Vec<Vec3<f64>> {
    let mut distance: Vec3<f64>;
    let mut acc: Vec<Vec3<f64>>;
    let mut radius: f64;
    let mut force: f64;

    // O(n^2)
    // Dude there is a lot of redundant code, binary sort this shit
    for (idx1, particle1) in particles.enumerate(){
        for (idx2, particle2) in particles[idx..].enumerate(){
            distance = particle1.pos - particle2.pos;
            radius = distance.norm();
            force =  field_strength * (particle1.mass + particle2.mass) / radius.powi(2);
            calc_accelerations[idx1] += -force / particle1.mass * distance / radius;
            calc_accelerations[idx2] += force / particle2.mass * distance / radius;
        }
    }
    return acc
}



/// Simulation builds to n body
/// 
/// 1 body - nada
/// 2 body - single for loop
/// 3 body - vectorized
/// n body - barnes hut tree

fn main() {

    // 1 body
    let mut particle1: Particle = Particle {
        pos = Vec3::new(0.,0.,0.),
        vel = Vec3::new(0.,0.,1.),
        mass = 3.0
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
    let mut distance: Vec3<f64>;
    let mut radius: f64;
    let mut force: f64;
    let mut dir: Vec3<f64>;
    let mut acc_1: Vec3<f64>;    
    let mut acc_2: Vec3<f64>;

    // Control duration and precision
    for dt in (0..t_end).iter(){

        distance = particle1.pos - particle2.pos;
        radius = distance.norm();
        force = field_strength * (particle1.mass + particle2.mass) / radius.powi(2);
        dir = distance / radius;

        acc_1 = -force / particle1.mass * dir;
        acc_2 = force / particle2.mass * dir;

        particle1.pos += particle1.vel * dt;
        particle2.pos += particle2.vel * dt;

        particle1.vel += acc1 * dt;
        particle2.vel += acc2 * dt;

    }

    // 3 body - initial vectorization
    let mut particle3: Particle = Particle {
        pos = Vec3::new(0.,0.,3.),
        vel = Vec3::new(2.,0.,0.),
        mass = 7.0
    };

    let mut particles: vec![particle1, particle2, particle3];

    for dt in (0..t_end).iter(){
        let accel: Vec<Vec3<f64>> = calc_acceleration(particles, field_strength);
        particles.zip(accel.iter()).for_each(|(particle, accel)| {
            particle.pos += particle.vel * dt;
            particle.vel += accel * dt;
        });

    }



    // N body- watch your Kolmogrov Complexity

}