/*
Barnes hut multibody propogator
*/

pub const TOLERANCE: f64 = 1e-16;

use nalgebra::Vector3;
use crate::donnager::{spacetime as xyzt, constants as cst};


#[derive(Clone, Debug, PartialEq)]
pub struct BhNode {
    pub child_base: usize,
    pub child_mask: u16,
    pub quad_size: f64,
    pub mass: f64
}

impl BhNode {

    /// Initialize new Node
    /// 
    /// Inputs
    /// ------
    /// range: `(Vector3<f64>, Vector3<f64>)`
    ///     Min and max position vectors
    /// 
    /// Outputs
    /// -------
    /// node: `Node`
    ///     Initialized `Node` struct
    fn new(range: (Vector3<f64>,Vector3<f64>)) -> Self {
        let node: BhNode = BhNode {
            child_base: 0,
            child_mask: 0,
            mass: 0.,
            quad_size: (range.0 - range.1).norm()
        };
        node
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BhTree {
    pub nodes: Vec<BhNode>,
    pub center_of_mass: Vec<Vector3<f64>>,
    pub avg_weighted_pos: Vec<Vector3<f64>>,
    pub mins: Vec<Vector3<f64>>,
    pub maxs: Vec<Vector3<f64>>,
    pub theta_sq: f64
}

impl BhTree {

    /// Initialize new `Tree`
    /// 
    /// Inputs
    /// ------
    /// particles: `Vector<Particle>`
    ///     List of particles to store
    /// 
    /// theta: `f64`
    ///     Angular precision
    /// 
    /// Outputs
    /// -------
    /// tree: `Tree`
    ///     Initialized `Tree` struct
    pub fn new(particles: Vec<xyzt::Particle>, theta: f64) -> BhTree {

        let range: (Vector3<f64>, Vector3<f64>) = calc_range(particles.clone());

        let mut tree: BhTree = BhTree {
            nodes: Vec::with_capacity(particles.len()),
            mins: Vec::with_capacity(particles.len()),
            maxs: Vec::with_capacity(particles.len()),
            center_of_mass: Vec::with_capacity(particles.len()),
			avg_weighted_pos: Vec::with_capacity(particles.len()),
            theta_sq: theta.powi(2)
        };

        tree.nodes.push(BhNode::new(range));
        tree.mins.push(range.0);
        tree.maxs.push(range.1);
		tree.avg_weighted_pos.push(Vector3::zeros());
        tree.add_particles_to_node(&particles, 0);

        for (idx, node) in &mut tree.nodes.iter().enumerate() {
			tree.center_of_mass.push(tree.avg_weighted_pos[idx] / node.mass);
		}

        return tree
    }


    /// Add list of Particles to node
    /// 
    /// Inputs
    /// ------
    /// particles: `Vector<&Particle>`
    ///     List of particles to store
    /// 
    /// node_id: `usize`
    ///     Node identification number. 
    pub fn add_particles_to_node(
        &mut self, 
        particles: &Vec<xyzt::Particle>, 
        node_id: usize
    ) {

        // If singular particle (end of tree branch)
        // node mass += particle mass 
        // av weight pos += particle av weight pos
        if particles.len() == 1 {
            self.nodes[node_id].mass += particles[0].mass;
			self.avg_weighted_pos[node_id] += particles[0].motion[0] * particles[0].mass;
            return;
        };

        // TODO-TD: improve initialization
        let mut particle_trees: [Vec<xyzt::Particle>; 8] = [
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
            Vec::with_capacity(particles.len() / 4),
        ];
        
        let center: Vector3<f64> = (self.mins[node_id] + self.maxs[node_id]) / 2.0;
        for particle in particles {
            self.nodes[node_id].mass += particle.mass;
			self.avg_weighted_pos[node_id] += particle.motion[0] * particle.mass;
            
            let offset: Vector3<f64> = particle.motion[0] - center;
            let x_offset: usize = if offset.x > 0.0 {0} else {1};
            let y_offset: usize = if offset.y > 0.0 {0} else {1};
            let z_offset: usize = if offset.z > 0.0 {0} else {1};
            let index: usize = x_offset + y_offset * 2 + z_offset * 4;

            particle_trees[index].push(particle.clone());
        };

        self.nodes[node_id].child_base = self.nodes.len();
        for (idx, particle_tree) in particle_trees.iter().enumerate() {
            if !particle_tree.is_empty() {
                let old_range: (Vector3<f64>, Vector3<f64>) = 
                    (self.mins[node_id], self.maxs[node_id]);
                let new_range: (Vector3<f64>, Vector3<f64>) = 
                    BhTree::get_bounding_box(old_range, center, idx);
                
                self.nodes.push(BhNode::new(new_range));
                self.mins.push(new_range.0);
                self.maxs.push(new_range.1);
				self.avg_weighted_pos.push(Vector3::zeros());
                self.nodes[node_id].child_mask |= 1 << idx;
            }
        }
        
        let mut child_node: usize = self.nodes[node_id].child_base;

        for particle_tree in particle_trees.iter() {
            if !particle_tree.is_empty() {
                self.add_particles_to_node(particle_tree, child_node);
                child_node += 1;
            }
        }
    }


    /// Get bounding box of Tree
    /// 
    /// Inputs
    /// ------
    /// range: `(Vector3<f64>, Vector3<f64>)`
    ///     Min and max points of space
    /// 
    /// center: `Vector3<f64>`
    ///     Center point of Tree
    /// 
    /// node_id: `usize`
    ///     Node identification number
    /// 
    /// Outputs
    /// -------
    /// range: `(Vector3<f64>, Vector3<f64>)`
    ///     Min and max position points
    pub fn get_bounding_box(
        range: (Vector3<f64>, Vector3<f64>), 
        center: Vector3<f64>, 
        node_id: usize
    ) -> (Vector3<f64>, Vector3<f64>) {
        let mut min_coord: Vector3<f64> = center;
		let mut max_coord: Vector3<f64> = range.1;

		if node_id == 1 {
			min_coord.x = range.0.x;
			max_coord.x = center.x;
		} 

		if node_id == 2 {
			min_coord.y = range.0.y;
			max_coord.y = center.y;
		} 

		if node_id == 4 {
			min_coord.z = range.0.z;
			max_coord.z = center.z;
		} 
		return (min_coord, max_coord)
	    }


    /// Calculate acceleration of a node
    /// 
    /// Inputs
    /// ------
    /// pos : `Vector3<f64>`
    ///     Position vector
    /// 
    /// node_id : `usize`
    ///     Node identification number
    /// 
    /// Outputs
    /// -------
    /// acc : `Vector3<f64>`
    ///     Acceleration vector
    pub fn barnes_hut_node_acc(
        &self, 
        pos: Vector3<f64>, 
        node_id: usize
    ) -> Vector3<f64> {
        let curr_node = &self.nodes[node_id];
        let dist: Vector3<f64> = self.center_of_mass[node_id] - pos;
        let dist_sq: f64 = dist.norm_squared();
        let is_close: bool = curr_node.quad_size < self.theta_sq * dist_sq;

        if is_close || curr_node.child_mask == 0 {
            if dist_sq < TOLERANCE {
                return Vector3::new(0.,0.,0.);
            };
            return dist * cst::GRAV_CONST * (curr_node.mass) / dist_sq;
        };

        let mut node: usize = curr_node.child_base;
        let mut sum: Vector3<f64> = Vector3::zeros();

        for node_idx in 0..8 {
            if curr_node.child_mask & (1 << node_idx) != 0 {
                sum += self.barnes_hut_node_acc(pos, node);
                node += 1;
            }
        }
        return sum;

        }

    }
 
/// Propogate set of gravitationally attracting particles
/// 
/// Inputs
/// ------
/// particles: `Vec<Particle>` 
///     Vector of Particles
/// 
/// step_size: `f64`
///     Time in seconds to increment simulation
/// 
/// n_steps: `usize`
///     Number of steps to propogate
/// 
/// theta: `f64`
///     Angular precision for barnes hut calculation
/// 
/// is_debug: `bool`
///     Show propogation
/// 
/// Outputs
/// -------
/// particles: `Vec<Particle>` 
///     Vector of Particles updated at t = step_size * n_steps 
pub fn barnes_hut_gravity(
    mut particles: Box<[xyzt::Particle]>,
    step_size: f64,
    n_steps: usize,
    theta: f64,
    is_debug: bool
) -> Box<[xyzt::Particle]> {
    let mut motion = vec![Vector3::zeros(); 3].into_boxed_slice();
    
    for i_step in 0..n_steps {
        // Create Tree
        let tree = BhTree::new(particles.to_vec(), theta);

        if is_debug{
            println!("t = {}", i_step as f64 * step_size)
        };

        // TODO-TD: multithread
        particles.iter_mut().for_each(|part| {
            motion[2] = tree.barnes_hut_node_acc(part.motion[0], 0);
            motion[1] = part.motion[1] + motion[2] * step_size;
            motion[0] = part.motion[0] + part.motion[1] * step_size + 0.5 * motion[2] * step_size * step_size;

            part.motion = motion.to_vec();
            
            if is_debug {
                print!("pos: {:.5?} \t\t", part.motion[0]);
            }
        });
    }
    particles
}

/// Calculate range of min and max positions of set of Particles
/// 
/// Inputs
/// ------
/// particles: `Vec<Particle>`
///     List of particles
/// 
/// Outputs
/// -------
/// range: `(Vector3<f64>, Vector3<f64>)`
///     Min and max points of space
pub fn calc_range(
    particles: Vec<xyzt::Particle>
) -> (Vector3<f64>, Vector3<f64>) {
    let mut min: Vector3<f64> = Vector3::zeros();
    let mut max: Vector3<f64> = Vector3::zeros();

    particles.iter()
        .for_each(|particle| {
            min.iter_mut()
                .zip(max.iter_mut())
                .zip(particle.motion[0].iter())
                .for_each(|((min_val, max_val), pos)| { 
                    if pos < min_val{*min_val=*pos;}
                    if pos > max_val{*max_val=*pos;}
                })
            }
        );
    return (min, max)
    }


#[cfg(test)]
mod barneshut_tests{
    use super::*;
    
    #[test]
    fn test_prop(){
        
        // Earth
        let earth: xyzt::Body = xyzt::Body {
            name: String::from("Earth"),
            grav_param: cst::EARTH::GRAV_PARAM,
            eq_radius: cst::EARTH::RADIUS_EQUATOR,
            rotation_rate: cst::EARTH::ROT_RATE,
            eccentricity: cst::EARTH::ECC
        };

        let radius: f64 = earth.calc_stationary_orbit();
        assert_eq!(radius, 42163779.55713436);

        let vel: f64 = earth.calc_orbital_velocity_mag(radius);
        assert_eq!(vel, 4348.18527478043);

        let radius_vec: Vector3<f64> = Vector3::new(radius, 0., 0.);
        let vel_vec: Vector3<f64> = Vector3::new(0., vel, 0.);
        let acc_vec: Vector3<f64> = earth.calc_body_grav(radius_vec);
        let sat_motion_0: Vec<Vector3<f64>> = vec![radius_vec, vel_vec, acc_vec];
        assert_eq!(acc_vec, Vector3::new(0.2242056497591453, 0., 0.));

        let pos_lla: Vector3<f64> = xyzt::ecef_to_lla(radius_vec);
        assert_eq!(pos_lla, Vector3::new(0.0, 0.0, 42157401420.13436));
        let satellite: xyzt::Particle = xyzt::Particle { mass: 5e4, motion: sat_motion_0};

        let motion: Vec<Vector3<f64>> = vec![Vector3::zeros(); 3];
        let earth_particle: xyzt::Particle = earth.to_particle(motion);
        assert_eq!(earth_particle.mass, 5.972e24);

        let mut particles: Box<[xyzt::Particle]> = vec![earth_particle, satellite].into_boxed_slice();
    
        let step_size: f64 = 0.1;
        let theta: f64 = 1.0;
        let n_steps: usize = 10;
        let is_debug: bool = false;
    
        particles = barnes_hut_gravity(particles, step_size, n_steps, theta, is_debug);
        let expected: Vector3<f64> = Vector3::new(37361232.85564361, 4201.215848944788, 0.0);
        assert_eq!(particles[1].motion[0],expected);
    }
}
