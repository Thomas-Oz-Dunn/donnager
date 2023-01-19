/*
Gravitational Bodies
*/

use nalgebra as na;
use std::f64::consts::PI;
use na::Vector3;

use crate::constants as cst;

pub const TOLERANCE: f64 = 1e-8;

#[derive(Clone, Debug, PartialEq)]
pub struct Body{
    pub name: String,
    pub grav_param: f64,
    pub eq_radius: f64,
    pub rotation_rate: f64,
    pub eccentricity: f64
}


impl Body {

    // Calculate gravitational acceleration at radial distance
    pub fn calc_grav_acc(&self, radius: f64) -> f64 {
        let grav_acc: f64 = self.grav_param / radius.powi(2);
        return grav_acc
    }
    

    // Calculate required orbital velocity at radial distance
    pub fn calc_orbital_velocity(&self, radius: f64) -> f64 {
        // TODO-TD: Vectorize
        let vel: f64 = (2.0 * self.grav_param / radius).sqrt();
        return vel
    }

    // Calculate period of orbit
    pub fn calc_period(&self, semi_major_axis: f64) -> f64 {
        let time: f64 = 2.0 * PI * (semi_major_axis.powi(3)/self.grav_param).sqrt();
        return time
    }
    
    // Calculate radius for stationary orbit above body surface
    pub fn calc_stationary_orbit(&self) -> f64 {
        let period: f64 = 2.0 * PI / self.rotation_rate;
        let a: f64 = self.grav_param * period.powi(2); // 
        let r_mag: f64 = (a / (4.0 * PI.powi(2))).powf(1.0 / 3.0); // 
        return r_mag
    }

    // Geodetic to rectangular coordinates
    // E.g. Latitude, Longitude, Altitude to ECEF
    pub fn geodetic_to_xyz(&self, lla: Vector3<f64>) -> Vector3<f64> {
        let radius: f64 = self.calc_prime_vertical(lla[0]);
        let x: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].cos();
        let y: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].sin();
        let z: f64 = ((1.0 - self.eccentricity.powi(2)) * radius + lla[2]) * lla[0].sin();
        let xyz: Vector3<f64> = Vector3::new(x, y, z); 
        return xyz
    }

    // Calculate prime vertical radius to surface at latitude
    pub fn calc_prime_vertical(&self, lat_deg: f64) -> f64 {
        let lat_radians: f64 = PI * lat_deg / 180.0;
        let radius: f64 = 
            self.eq_radius / (1.0 - (self.eccentricity * lat_radians.sin()).powi(2)).sqrt();
        return radius
    }

    // Rectangular coordinates to geodetic
    // E.g. ECEF to LLH
    pub fn xyz_to_geodetic(&self, xyz: Vector3<f64>) -> Vector3<f64> {
        // Zhu's method
        let a: f64 = self.eq_radius;
        let ecc_2: f64 = self.eccentricity.powi(2);

        let b: f64 = (a.powi(2)*(1.0 - ecc_2)).sqrt();
        let ecc_2_prime: f64 = a.powi(2) / b.powi(2) - 1.0;
        let p: f64 = (xyz[0].powi(2) + xyz[1].powi(2)).sqrt();
        let g: f64 = p.powi(2) + (1.0 - ecc_2) * xyz[2].powi(2) - 
            ecc_2 * (a.powi(2) - b.powi(2));
        let c: f64 = ecc_2.powi(2) * 54.0 * b.powi(2) * xyz[2].powi(2) * p.powi(2) / (g.powi(3));
        let s: f64 = (1.0 + c + (c.powi(2) + 2.0 * c).sqrt()).powf(1.0 / 3.0);
        let capP: f64 = 54.0 * b.powi(2)*xyz[2].powi(2) / 
            (3.0 * (s + 1.0 + 1.0 / s).powi(2) * g.powi(2));
        let q: f64 = (1.0 + 2.0 * ecc_2.powi(2) * capP).sqrt();
        let r_0: f64 = -capP * ecc_2 * p /(1.0+q) + 
            ((a.powi(2)/2.0)*(1.0 + 1.0 / q) - 
            capP * (1.0 - ecc_2) * xyz[2].powi(2) / (q * (1.0 + q)) - 
            capP*p.powi(2)/2.0).sqrt();
        let u: f64 = ((p - ecc_2*r_0).powi(2) + xyz[2].powi(2)).sqrt();
        let v: f64 = ((p - ecc_2*r_0).powi(2) + (1.0 - ecc_2)*xyz[2].powi(2)).sqrt();
        let z_0: f64 = b.powi(2) * xyz[2] / (a * v);

        let alt: f64 = u * (1.0 - b.powi(2) / (a * v));
        let lat: f64 = ((xyz[2] + ecc_2_prime*z_0)/p).atan();
        let lon: f64 = (xyz[1] / xyz[0]).atan();
        let lla: Vector3<f64> = Vector3::new(lat, lon, alt);
        return lla
    }

}

#[derive(Clone, Debug, PartialEq)]
pub struct Particle {
    pub mass: f64,
    pub motion: Vec<Vector3<f64>>,
}


impl Particle {
   
    pub fn update(
        mut self, 
        mut particle2: Particle,
        time_step: f64
    ) -> Particle {
        let mut delta : Vector3<f64>; 
        delta = (self.motion[1] + 0.5 * self.motion[2] * time_step) * time_step;
        self.motion[0] += delta;
        
        delta = (particle2.motion[1] + 0.5 * particle2.motion[2] * time_step) * time_step;
        particle2.motion[0] += delta;
        
        delta = particle2.motion[2] * time_step;
        self.motion[1] += delta;

        delta = particle2.motion[2] * time_step;
        particle2.motion[1] += delta;


        let distance: Vector3<f64> = self.motion[0] - particle2.motion[0];
        let radius: f64 = distance.norm();
        let force: f64 = cst::GRAV_CONST * (self.mass + particle2.mass) / radius.powi(2);
        let dir: Vector3<f64> = distance / radius;
        
        self.motion[2] += -force / self.mass * dir;
        particle2.motion[2] += force / particle2.mass * dir;
        particle2
    }

}


/// Calculate acceleration, brute force method (O(n^2) runtime)
/// 
/// Inputs
/// ------
/// particles: `Vec<Particle>` [N, ]
///     Vector of `Particle` structs to be simulated
/// 
/// field_strength: `f64`
///     Scalar value describing field of interest. 
///     Gravity:
///     Electromagnetism:
/// 
/// Outputs
/// -------
/// accelerations: `Vec<Vector3<f64>>` [N, 3]
///     Net accleration vector for each particle
pub fn calc_acceleration(
    particles: Vec<Particle>
) -> Vec<Vector3<f64>> {

    let mut distance: Vector3<f64>;
    let mut radius: f64;
    let mut force: f64;
    let zero_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    let mut acc_mut: Vec<Vector3<f64>> = vec![zero_vec; particles.len()];

    // O(n^2)
    // Dude there is a lot of redundant code, binary sort this shit
    for (idx1, particle1) in particles.iter().enumerate(){
        for (idx2, particle2) in particles[idx1..particles.len()].iter().enumerate(){
            distance = particle1.motion[0] - particle2.motion[0];
            radius = distance.norm();
            force =  cst::GRAV_CONST * (particle1.mass + particle2.mass) / radius.powi(2);
            acc_mut[idx1] -= force / particle1.mass * distance / radius;
            acc_mut[idx2] += force / particle2.mass * distance / radius;
        }
    }
    let acc: Vec<Vector3<f64>> = acc_mut;
    return acc
}



pub struct Node {
    pub child_base: usize,
    pub child_mask: u16,
    pub quad_size: f64,
    pub mass: f64
}

impl Node {
    fn new(min: Vector3<f64>, max: Vector3<f64>) -> Self {
        let q_size = (max - min).abs();
        Node {
            child_base: 0,
            child_mask: 0,
            mass: 0,
            quad_size: q_size
        }
    }
}

pub struct Tree {
    pub nodes: Vec<Node>,
    pub center_of_mass: Vec<Vector3<f64>>,
    pub mins: Vec<Vector3<f64>>,
    pub maxs: Vec<Vector3<f64>>,
    pub theta_sq: f64
}

impl Tree {
    // Tree from vector of particles

    pub fn new(particles: Vec<Particle>, theta: f64) -> Self {
        let range = find_range(particles);
        let mut tree = Tree::empty(range, theta, particles.len());
        tree.add_particles_to_node(particles, 0);
        tree
    }

    pub fn add_particles_to_node(
        self, 
        particles: Vec<Particle>, 
        node: usize
    ) {
        if particles.len() == 1 {
            self.nodes[node].mass += particles[0].mass;
            return self;
        };

        let center_point = (self.mins[node] + self.maxs[node]) / 2.0;

        let mut particle_trees: [Vec<&Particle>; 8] = [
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
        ];

        for &particle in particles {
            self.nodes[node].mass += particle.mass;
            let index = Tree::get_id_from_center(center_point, particle.position);
            particle_trees[index].push(particle);
        }

        //Recurse
        self.nodes[node].child_base = self.nodes.len();
        for (idx, particle_tree) in particle_trees.iter().enumerate() {
            if !particle_tree.is_empty() {
                let (min, max) = Tree::get_bounding_box_from_id(
                    self.mins[node],
                    self.maxs[node],
                    center_point,
                    idx,
                );

                self.nodes.push(Node::new(min, max));
                self.mins.push(min);
                self.maxs.push(max);
                self.nodes[node].children_mask |= 1 << i;
            }
        }

        let mut child_node = self.nodes[node].child_base;
        for particle_tree in particle_trees.iter() {
            if !particle_tree.is_empty() {
                self.add_particles_to_node(particle_tree, child_node);
                child_node += 1;
            }
        }
    }

    pub fn empty(
        range: (Vector3<f64>,Vector3<f64>), 
        theta: f64, 
        n_nodes: usize
    ) -> Self {
        let mut tr = Tree {
            nodes: Vec::with_capacity(n_nodes),
            mins: Vec::with_capacity(n_nodes),
            maxs: Vec::with_capacity(n_nodes),
            center_of_mass: Vec::with_capacity(n_nodes),
            theta_sq: theta.powi(2)
        };
        tr.nodes.push(Node::new(range.0, range.1));
        tr.mins.push(range.0);
        tr.maxs.push(range.1);
        return tr
    }

    pub fn get_bounding_box_from_id(
        min: Vector3<f64>, 
        max: Vector3<f64>, 
        center: Vector3<f64>, 
        idx: usize
    ) -> (Vector3<f64>, Vector3<f64>) {

    }


}


pub fn barnes_hut_gravity(
    particles: Vec<Particle>
){

}


pub fn find_range(particles: Vec<Particle>) -> (Vector3<f64>, Vector3<f64>) {
    let mut min: Vector3<f64> = Vector3::zeros();
    let mut max: Vector3<f64> = Vector3::zeros();

    particles.iter().for_each(
        |particle| {
            particle.motion[0].iter_mut().zip(
                min.iter_mut().zip(max.iter_mut())).for_each(
                |(pos, (min_val, max_val))| { 
                    if pos < min_val{min_val=pos;}
                    if pos < max_val{max_val=pos;}
                })
            }
        );
    return (min, max)
    }
