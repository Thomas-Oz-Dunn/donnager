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
        mass: 30.0,
        motion: vec![
            Vector3::new(0.,2.8,0.),
            Vector3::zeros(),
            Vector3::zeros()]
    };

    // Body
    let particle2: grav::Particle = grav::Particle {
        mass: 50.0,
        motion: vec![
            Vector3::new(0.7,0.5,0.),
            Vector3::zeros(),
            Vector3::zeros()]
    };

    // Body 
    let particle3: grav::Particle = grav::Particle {
        mass: 15.0,
        motion: vec![
            Vector3::new(0.,0.,1.),
            Vector3::zeros(),
            Vector3::zeros()]
    };
    
    // N body- watch your Kolmogrov Complexity
    let theta: f64  = 1.0;
    let n_steps: usize = 1000;
    let step_size: f64 = 0.01;
    let is_debug: bool = true;

    let mut particles: Vec<grav::Particle> = [particle1, particle2, particle3].to_vec();

    println!("\nStart\n");
    for (particle, idx) in particles.iter().zip(0..){
        if is_debug {
            print!("{} pos: {:.5} x, {:.5} y, {:.5} z\t\t", idx, 
                particle.motion[0][0], particle.motion[0][1], particle.motion[0][2]);
        }
    }

    particles = grav::barnes_hut_gravity(
        particles, 
        step_size, 
        n_steps, 
        theta, 
        is_debug);

}