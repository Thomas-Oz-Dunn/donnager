/*
n-body simulation demonstration
Barnes Hut
*/

use std::fs;
use std::ops::Range;

use plotters::prelude::*;
use plotters::coord::types::*;
use nalgebra as na;
use na::Vector3;


#[derive(Clone)]
pub struct Particle {
    pub pos: Vector3<f64>,
    pub vel: Vector3<f64>,
    pub mass: f64,
}

pub struct Node {
    pub mass: f64
}

impl Node {

}

pub struct Tree {
    pub nodes: Vec<Node>,
    pub center_of_mass: Vector3<f64>
}

impl Tree {

}


/// Calculate acceleration, brute force
fn calc_acceleration(
    particles: Vec<Particle>,
    field_strength: f64
) -> Vec<Vector3<f64>> {

    let mut distance: Vector3<f64>;
    let mut radius: f64;
    let mut force: f64;
    let mut acc_mut: Vec<Vector3<f64>> = Vec::with_capacity(particles.len());

    // O(n^2)
    // Dude there is a lot of redundant code, binary sort this shit
    for (idx1, particle1) in particles.iter().enumerate(){
        for (idx2, particle2) in particles[idx1..].iter().enumerate(){
            distance = particle1.pos - particle2.pos;
            radius = distance.norm();
            force =  field_strength * (particle1.mass + particle2.mass) / radius.powi(2);

            acc_mut[idx1] += -force / particle1.mass * distance / radius;
            acc_mut[idx2] += force / particle2.mass * distance / radius;
        }
    }
    let acc: Vec<Vector3<f64>> = acc_mut;
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
        pos: Vector3::new(0.,0.,0.),
        vel: Vector3::new(0.,0.,1.),
        mass: 3.0
    };

    // 2 body
    let mut particle2: Particle = Particle {
        pos: Vector3::new(1.,0.,0.),
        vel: Vector3::new(0.,-2.,0.),
        mass: 5.0
    };

    match fs::create_dir("./images") {
        Err(why) => println!("! {:?}", why.kind()),
        Ok(_) => {},
    }
    let root_drawing_area = BitMapBackend::gif(
        "./images/animated.gif", 
        (500, 500), 
        1_000  /* Each frame show 1s */
    ).unwrap().into_drawing_area();

    root_drawing_area.fill(&WHITE).unwrap();
    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, 30)
        .set_label_area_size(LabelAreaPosition::Bottom, 30)
        .build_cartesian_2d::<Range<i32>, Range<i32>>(0..100, 0..100)
        .unwrap();

    ctx.configure_mesh().draw().unwrap();

    // Run
    let t_end: i32 = 60;
    let field_strength: f64 = 4.0;
    let mut distance: Vector3<f64>;
    let mut radius: f64;
    let mut force: f64;
    let mut dir: Vector3<f64>;
    let mut acc_1: Vector3<f64>;    
    let mut acc_2: Vector3<f64>;

    // Control duration and precision
    for _ in [0..t_end].iter(){

        distance = particle1.pos - particle2.pos;
        radius = distance.norm();
        force = field_strength * (particle1.mass + particle2.mass) / radius.powi(2);
        dir = distance / radius;
        
        acc_1 = -force / particle1.mass * dir;
        acc_2 = force / particle2.mass * dir;
        
        particle1.pos += particle1.vel;
        particle2.pos += particle2.vel;

        particle1.vel += acc_2;
        particle2.vel += acc_2;

    }

    // 3 body - initial vectorization
    let mut particle3: Particle = Particle {
        pos: Vector3::new(0.,0.,3.),
        vel: Vector3::new(2.,0.,0.),
        mass: 7.0
    };

    let mut particles: Vec<Particle> = Vec::new();
    particles = [particle1, particle2, particle3].to_vec();

    for _ in [0..t_end].iter(){
        
        let accel: Vec<Vector3<f64>> = calc_acceleration(particles.clone(), field_strength);
        (0..).zip(
            accel.iter().zip(
                particles.iter()))
                .map(|(i_point, (acc, particle))| {
            
                ctx.draw_series(LineSeries::new(
                    [particle.pos[0], particle.pos[1]].iter(), &Palette99::pick(i_point)
                )).unwrap();
            
                particle.pos += particle.vel;
                particle.vel += acc;
            });

    }



    // N body- watch your Kolmogrov Complexity

}