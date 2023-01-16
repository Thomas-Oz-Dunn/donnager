/*
n-body simulation demonstration
Barnes Hut
*/

use std::fs;
use std::ops::Range;
use plotters::prelude::*;
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

    // Render plot
    match fs::create_dir("./images") {
        Err(why) => println!("! {:?}", why.kind()),
        Ok(_) => {},
    }

    let root_drawing_area = BitMapBackend::gif(
        "./images/animated.gif", 
        (500, 500), 
        1_000  /* Each frame show 1s */
    ).unwrap().into_drawing_area();

    let x_spec: Range<f64> = -10.0..10.0;
    let y_spec: Range<f64> = -10.0..10.0;

    root_drawing_area.fill(&WHITE).unwrap();
    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, 30)
        .set_label_area_size(LabelAreaPosition::Bottom, 30)
        .build_cartesian_2d::<Range<f64>, Range<f64>>(x_spec, y_spec)
        .unwrap();

    ctx.configure_mesh().draw().unwrap();


    // Run
    // 1 body
    let particle1: grav::Particle = grav::Particle {
        mass: 30.0,
        pos: Vector3::new(0.9,2.8,0.),
        vel: Vector3::new(0.2,-0.1, 0.),
        acc: Vector3::zeros()
    };

    // 2 body
    let mut particle2: grav::Particle = grav::Particle {
        mass: 50.0,
        pos: Vector3::new(0.7,0.5,0.),
        vel: Vector3::new(-0.1,0.07,0.),
        acc: Vector3::zeros()
    };

    let t_end: i32 = 60;

    // Control duration and precision
    for _ in [0..t_end].iter(){
        particle2 = particle1.update(particle2, 1.0);
    }



    // 3 body - initial vectorization
    let particle3: grav::Particle = grav::Particle {
        mass: 0.0,
        pos: Vector3::new(0.,0.,0.),
        vel: Vector3::new(2.,0.,0.),
        acc: Vector3::new(0.,0., 0.),
    };

    let mut particles: Vec<grav::Particle> = [particle1, particle2, particle3].to_vec();

    for _ in [0..t_end].iter(){
        
        let accel: Vec<Vector3<f64>> = grav::calc_acceleration(particles.clone());
        (0..).zip(
            accel.iter().zip(
                particles.iter_mut()))
                    .for_each(|(i_point, (acc, particle))| {
                        ctx.draw_series(
                            // PLot idea, circle at point!
                            LineSeries::new(
                                [(particle.pos[0], particle.pos[1]), (particle.pos[0] + particle.vel[0], particle.pos[1] + particle.vel[1])], 
                                Palette99::pick(i_point))
                            ).unwrap().label(format!("Particle {}", i_point));

                particle.pos += particle.vel;
                particle.vel += acc;
            });

    }
    ctx.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();



    // N body- watch your Kolmogrov Complexity

}