/*
n-body simulation demonstration using Barnes Hut method
*/

use nalgebra as na;
use na::Vector3;

use donnager::gravity as grav;
use donnager::cosmos as cosm;

fn main() {

    // Body
    let particle1: cosm::space::Particle = cosm::space::Particle {
        mass: 1e13,
        motion: vec![
            Vector3::new(5.,1.,0.),
            Vector3::zeros(),
            Vector3::zeros()]
    };

    // Body
    let particle2: cosm::space::Particle = cosm::space::Particle {
        mass: 1e13,
        motion: vec![
            Vector3::new(7.,4.,0.),
            Vector3::zeros(),
            Vector3::zeros()]
    };

    // Body 
    let particle3: cosm::space::Particle = cosm::space::Particle {
        mass: 1e13,
        motion: vec![
            Vector3::new(0.,10.,0.),
            Vector3::zeros(),
            Vector3::zeros()]
    };
    
    // N body configuration
    let theta: f64  = 1.0;
    let n_steps: usize = 66;
    let step_size: f64 = 0.01;
    let is_debug: bool = true;

    let mut particles  = vec![particle1, particle2, particle3].into_boxed_slice();

    println!("\nStart");

    for particle in particles.iter(){
        print!("pos: {:.5?} \t\t", particle.motion[0]);
    }

    particles = grav::barneshut::barnes_hut_gravity(
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
