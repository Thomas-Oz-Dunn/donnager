/*
n-body simulation demonstration
Barnes Hut
*/

use nalgebra as na;
use na::Vector3;

use donnager::cosmos::grav as grav;

fn main() {

    // Body
    let particle1: grav::Particle = grav::Particle {
        mass: 1e13,
        motion: vec![
            Vector3::new(5.,1.,0.),
            Vector3::zeros(),
            Vector3::zeros()]
    };

    // Body
    let particle2: grav::Particle = grav::Particle {
        mass: 1e13,
        motion: vec![
            Vector3::new(7.,4.,0.),
            Vector3::zeros(),
            Vector3::zeros()]
    };

    // Body 
    let particle3: grav::Particle = grav::Particle {
        mass: 1e13,
        motion: vec![
            Vector3::new(0.,10.,0.),
            Vector3::zeros(),
            Vector3::zeros()]
    };
    
    // N body configuration
    let theta: f64  = 1.0;
    let n_steps: usize = 100;
    let step_size: f64 = 0.01;
    let is_debug: bool = false;

    let mut particles: Vec<grav::Particle> = [particle1, particle2, particle3].to_vec();

    println!("\nStart");

    for particle in particles.iter(){
        print!("pos: {:.5?} \t\t", particle.motion[0]);
    }

    particles = grav::barnes_hut_gravity(
        particles, 
        step_size, 
        n_steps, 
        theta, 
        is_debug);

        
    println!("\nEnd @ t = {} s", (n_steps as f64) * step_size);
    for particle in particles.iter(){
        print!("pos: {:.5?} \t\t", particle.motion[0]);
    }

}