/*
n-body simulation demonstration
Barnes Hut
*/

use std::ops::Range;
use nalgebra as na;
use na::Vector3;

use donnager::cosmos::grav as grav;

/// Simulation builds to n body
/// 
/// 1 body - nada
/// 2 body - single for loop
/// 3 body - vectorized
/// n body - barnes hut tree

fn main() {

    // Run
    // 1 body
    let particle1: grav::Particle = grav::Particle {
        mass: 30.0,
        motion: vec![Vector3::new(0.9,2.8,0.),
                     Vector3::new(0.2,-0.1, 0.),
                     Vector3::zeros()]
    };

    // 2 body
    let mut particle2: grav::Particle = grav::Particle {
        mass: 50.0,
        motion: vec![Vector3::new(0.7,0.5,0.),
                     Vector3::new(-0.1,0.07,0.),
                     Vector3::zeros()]
    };

    // 3 body 
    let particle3: grav::Particle = grav::Particle {
        mass: 0.0,
        motion: vec![Vector3::new(0.,0.,0.),
                     Vector3::new(2.,0.,0.),
                     Vector3::new(0.,0., 0.)]
    };

    // N body- watch your Kolmogrov Complexity
    let mut particles: Vec<grav::Particle> = [particle1, particle2, particle3].to_vec();
    println!("\nStart\n");
    for (particle, idx) in particles.iter().zip(0..){
        println!("Particle {} pos: {} x, {} y, {} z \t", idx, 
            particle.motion[0][0], particle.motion[0][1], particle.motion[0][2])
    }
    let theta: f64  = 0.1;
    let n_steps: usize = 100;
    let step_size: f64 = 2.5;
    let is_debug: bool = true;

    particles = grav::barnes_hut_gravity(
        particles, 
        step_size, 
        n_steps, 
        theta, 
        is_debug);

    println!("\nEnd\n");
    for (particle, idx) in particles.iter().zip(0..){
        println!("\nParticle {} pos: {} x, {} y, {} z \n", idx, 
            particle.motion[0][0], particle.motion[0][1], particle.motion[0][2])
    }

}