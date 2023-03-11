/*
Gravitational Bodies
*/

use nalgebra as na;
use na::{Vector3, Matrix3};
use chrono::{DateTime, NaiveDateTime, NaiveDate, NaiveTime, TimeZone, Utc};
use std::f64::consts::PI;

use crate::{cosmos::spacetime as spacetime, constants as cst};

/// Orbit structure
/// 
/// Variables
/// ---------
/// name : `String`
#[derive(Clone, Debug, PartialEq)]
pub struct Orbit{
    pub name: String,
    pub grav_param: f64,
    pub semi_major_axis: f64, 
    pub eccentricity: f64,
    pub inclination: f64,
    pub raan: f64,
    pub argument_of_perigee: f64,
    pub mean_anomaly: f64,
    pub mean_motion: f64,
    pub epoch: DateTime<Utc>
}

impl Orbit {
    
    /// Populate Orbit from Keplerian parameters
    /// 
    /// Inputs
    /// ------
    /// name : `str`
    ///     Name of body this orbit is for
    /// 
    /// grav_param : `f64`           
    ///     Gravitational parameter of central body
    /// 
    /// semi_major_axis : `f64`
    ///     Semi-major axis of orbit in meters
    /// 
    /// eccentricity : `f64`  
    ///     Eccentricity of orbit, 0 <= eccentricity < 1.0
    /// 
    /// inclination : `f64`
    ///     Inclination of orbit in radians, 0 <= inclination < pi/2.0
    /// 
    /// raan : `f64`
    ///     Right ascension of the ascending node in radians, 0 <= raan < 2pi.
    /// 
    /// argument_of_perigee : `f64`
    ///     Argument of perigee in radians, 0 <= argument_of_perigee < 2pi.
    /// 
    /// mean_anomaly : `f64`
    ///     Mean anomaly in radians, 0 <= mean_anomaly < 2pi.
    /// 
    /// mean_motion : `f64`
    ///     Mean motion in radians per second.
    /// 
    /// epoch : `DateTime<Utc>`
    ///     Epoch of the orbit in UTC.
    /// 
    /// Outputs       
    /// -------
    /// orbit : `Orbit`           
    ///     Orbit structure with populated Keplerian parameters.
    pub fn from_keplerian(
        name: String,
        grav_param: f64,
        semi_major_axis: f64, 
        eccentricity: f64,
        inclination: f64,
        raan: f64,
        argument_of_perigee: f64,
        mean_anomaly: f64,
        mean_motion: f64,
        epoch: DateTime<Utc>
    ) -> Self {
        Orbit {
            name,
            grav_param,
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            argument_of_perigee,
            mean_anomaly,
            mean_motion,
            epoch
        }
    }

    /// Populate Orbit from standard Two Line Element
    /// 
    /// Inputs
    /// ------
    /// tle_str : `String` 
    ///     NORAD Two Line Element Identification String
    pub fn from_tle(
        tle_str: String
    ) -> Self {
        let lines: Vec<&str> = tle_str.lines().collect();
        let name: &str = lines[0];
        let bind1 = lines[1].to_string();
        let line1: Vec<&str> = bind1
            .split_whitespace()
            .collect();

        let epoch_str: &str = line1[3];
        let epoch_year: i32 = epoch_str[..=1]
            .to_string()
            .parse::<i32>()
            .unwrap();

        let year: i32;
        if epoch_year < 57{
            year = 2000 + epoch_year;
        } else {
            year = 1900 + epoch_year;
        }

        let binding = epoch_str[2..]
            .to_string();
        let epoch_day_full: Vec<&str> = binding
            .split_terminator('.')
            .collect();

        let day_of_year: u32 = epoch_day_full[0]
            .to_string()
            .parse::<u32>()
            .unwrap();

        let md: (u32, u32) = spacetime::calc_month_day(day_of_year, year);
        let date: NaiveDate = NaiveDate::from_ymd_opt(
            year, md.0, md.1).unwrap();
        
        let percent_of_day: f64 = 
        (".".to_owned() + &epoch_day_full[1].to_string())
            .parse::<f64>()
            .unwrap();

        let hours_dec: f64 = percent_of_day * 24.0;
        let hours_whole: u32 = hours_dec.div_euclid(24.0).floor() as u32;
        let hours_part: f64 = hours_dec.rem_euclid(24.0);
        
        let minutes_dec: f64 = hours_part * 60.;
        let minutes_whole: u32 = minutes_dec.div_euclid(60.).floor() as u32;
        let minutes_part: f64 = minutes_dec.rem_euclid(60.);

        let seconds_dec: f64 = minutes_part * 60.;
        let seconds_whole: u32 = seconds_dec.div_euclid(60.).floor() as u32;

        let time: NaiveTime = NaiveTime::from_hms_opt(
            hours_whole, minutes_whole, seconds_whole).unwrap();
        let dt: NaiveDateTime = NaiveDateTime::new(date, time);
        let epoch_date_time: DateTime::<Utc> = DateTime::<Utc>::from_utc(dt, Utc); 
        
        // let mean_motion_prime: &str = line1[4];
        // let mean_motion_2: &str = line1[5];
        
        let binding: String = lines[2].to_string();
        let line2: Vec<&str> = binding.split_whitespace().collect();
        
        let inc: f64 = line2[2]
            .to_string()
            .parse::<f64>()
            .unwrap();

        let raan: f64 = line2[3]
            .to_string()
            .parse::<f64>()
            .unwrap();

        let ecc: f64 =
        (".".to_owned() + & line2[4].to_string())
            .parse::<f64>()
            .unwrap();

        let arg_perigee: f64 = line2[5]
            .to_string()
            .parse::<f64>()
            .unwrap();

        let mean_anomaly: f64 = line2[6]
            .to_string()
            .parse::<f64>()
            .unwrap();

        let end_str: &str = line2[line2.len()-1];
        let mean_motion: f64 = end_str[..11]
            .to_string()
            .parse::<f64>()
            .unwrap();

        // Two Line element usage assumes Earth Centered
        let semi_major_axis: f64 = calc_semi_major_axis(
            cst::EARTH_GRAV_PARAM, mean_motion);
    
        Orbit {
            name: name.to_string(),
            grav_param: cst::EARTH_GRAV_PARAM,
            semi_major_axis,
            raan,
            eccentricity: ecc,
            inclination: inc,
            argument_of_perigee: arg_perigee,
            mean_anomaly,
            mean_motion,
            epoch: epoch_date_time
        }
    
    }

    /// Populate Orbit from position and velocity vectors
    /// 
    /// Inputs
    /// ------
    /// name : `str`
    ///     Name of object
    /// 
    /// grav_param : `f64`
    ///     Central body gravitational parameter
    /// 
    /// pos : `Vector3<f64>`
    ///     Position vector of object
    /// 
    /// vel : `Vector3<f64>`
    ///     Velocity vector of object
    /// 
    /// epoch_datetime: `DateTime<Utc>`
    ///     Epoch of object position and velocity vectors
    /// 
    /// Outputs	
    /// ------
    /// orbit : `Orbit`
    ///     Orbit object with populated fields.
    pub fn from_pos_vel(
        name: String,
        grav_param: f64,
        pos: Vector3<f64>,
        vel: Vector3<f64>,
        epoch_datetime: DateTime<Utc>
    ) -> Self {
        let spec_ang_moment: Vector3<f64> = pos.cross(&vel);
        let ecc_vec: Vector3<f64> = calc_ecc_vec(pos, vel, grav_param);
        let ascend_node_vec: Vector3<f64> = Vector3::z_axis().cross(&spec_ang_moment);

        let semi_major_axis: f64 = 
            spec_ang_moment.norm_squared() * (1.0 - ecc_vec.norm_squared()) / grav_param;

        Orbit {
            name,
            grav_param,
            semi_major_axis,
            eccentricity: ecc_vec.norm(),
            raan: calc_raan(ascend_node_vec),
            inclination: calc_inclination(spec_ang_moment),
            argument_of_perigee: ascend_node_vec.angle(&ecc_vec),
            mean_anomaly: ecc_vec.angle(&pos),
            mean_motion: calc_mean_motion(semi_major_axis, grav_param),
            epoch: epoch_datetime
        }

    }

    /// Calculate position and velocity vectors in perfocal frame
    /// 
    /// Inputs	
    /// ------
    /// time: f64
    ///     
    /// frame: `str`
    ///     Reference frame
    pub fn calc_pos_vel(
        &self, 
        time: f64, 
        frame: spacetime::ReferenceFrames
    ) -> (Vector3<f64>, Vector3<f64>) {
        let ecc: f64 = self.eccentricity;

        let mean_anom: f64 = (
            self.mean_anomaly + self.mean_motion * time) * cst::DEG_TO_RAD;
        let cos_mean_anom: f64 = mean_anom.cos();

        let ecc_anom: f64 = (
            mean_anom - ecc * cst::DEG_TO_RAD * (1.0 - cos_mean_anom)) * cst::DEG_TO_RAD;
        
        let true_anomaly_rad: f64 = 2.0 * (
            ecc_anom.sin() - ecc_anom.cos()).atan();
        let cos_true_anom: f64 = true_anomaly_rad.cos();
        let sin_true_anom: f64 = true_anomaly_rad.sin();

        let radius: f64 = self.semi_major_axis * (1.0 - ecc.powi(2)) / (1.0 + ecc * cos_true_anom);
        
        // Perifocal
        let x_pos: f64 = radius * cos_true_anom;
        let y_pos: f64 = radius * sin_true_anom;
        let z_pos: f64 = 0.0;
        let pos: Vector3<f64> = Vector3::new(x_pos, y_pos, z_pos);

        // Perifocal
        let x_vel: f64 = -self.mean_motion * radius * sin_true_anom;
        let y_vel: f64 = self.mean_motion * radius * (ecc + cos_true_anom);
        let z_vel: f64 = 0.0;
        let vel: Vector3<f64> = Vector3::new(x_vel, y_vel, z_vel);

        match frame {
            spacetime::ReferenceFrames::ECI => {
                let rotam: Matrix3<f64> = self.calc_pfcl_eci_rotam();
                let eci_pos: Vector3<f64> = rotam * pos;
                let eci_vel: Vector3<f64> = rotam * vel;
                return (eci_pos, eci_vel);
            },
            spacetime::ReferenceFrames::ECEF => {
                let pfcl_eci_rotam: Matrix3<f64> = self.calc_pfcl_eci_rotam();
                let eci_pos: Vector3<f64> = pfcl_eci_rotam * pos;
                let eci_vel: Vector3<f64> = pfcl_eci_rotam * vel;

                let new_time: f64 = self.epoch.timestamp() as f64 + time;
                let new_epoch_datetime: DateTime<Utc> = Utc.timestamp_opt(
                    new_time as i64, 0).unwrap();
                let eci_ecef_rotam: Matrix3<f64>  = spacetime::calc_eci_ecef_rotam(new_epoch_datetime);

                let ecef_pos: Vector3<f64> = eci_ecef_rotam * eci_pos;
                let ecef_vel: Vector3<f64> = eci_ecef_rotam * eci_vel;
                return (ecef_pos, ecef_vel);
            },
            spacetime::ReferenceFrames::PFCL => {
                return (pos, vel);
            }
        }
    }

    /// Calculate perifocal to eci rotation matrix
    pub fn calc_pfcl_eci_rotam(&self) -> Matrix3<f64> {
        let cos_raan: f64 = self.raan.cos();
        let sin_raan: f64 = self.raan.sin();
        let cos_inc: f64 = self.inclination.cos();
        let sin_inc: f64 = self.inclination.sin();
        let cos_arg_peri: f64 = self.argument_of_perigee.cos();
        let sin_arg_peri: f64 = self.argument_of_perigee.sin();

        let rot_mat: Matrix3<f64> = 
            Matrix3::new(cos_raan, -sin_raan, 0.0,
                             sin_raan, cos_raan, 0.0,
                             0.0, 0.0, 1.0);

        let rot_mat_2: Matrix3<f64> =
            Matrix3::new(cos_inc, 0.0, sin_inc,
                 0.0, 1.0, 0.0, 
                 -sin_inc, 0.0, cos_inc);

        let rot_mat_3: Matrix3<f64> = 
            Matrix3::new(cos_arg_peri, sin_arg_peri, 0.0, 
                -sin_arg_peri, cos_arg_peri, 0.0, 
                0.0, 0.0, 1.0);

        return rot_mat_3 * rot_mat_2 * rot_mat; 
    }

    /// Propogate orbit an increment of time
    /// 
    /// Inputs
    /// ------
    /// dt: `f64`
    ///     Time increment in seconds           
    /// 
    /// Outputs
    /// -------
    /// orbit: `Orbit`
    ///     Propogated orbit struct
    pub fn propogate(&self, dt: f64) -> Orbit {
        let mut new_orbit = self.clone();
        new_orbit.propogate_in_place(dt);
        new_orbit
    }
       
    /// Propogate orbit in place, without returning a new orbit instance.
    /// 
    /// Inputs
    /// ------
    /// dt: `f64`
    ///     Time step, in seconds.
    /// 
    /// Outputs
    /// -------
    /// None.
    pub fn propogate_in_place(&mut self, dt: f64) {
        let new_time: f64 = self.epoch.timestamp() as f64 + dt;
        let frame = spacetime::ReferenceFrames::ECI;
        let motion: (Vector3<f64>, Vector3<f64>) = self.calc_pos_vel(new_time, frame);
        let new_epoch_datetime: DateTime<Utc> = Utc.timestamp_opt(
            new_time as i64, 0).unwrap();

        let new_orbit: Orbit = Orbit::from_pos_vel(
            self.name.clone(), 
            self.grav_param, 
            motion.0, 
            motion.1, 
            new_epoch_datetime);
        *self = new_orbit;

        }
    }

/// Calculate the eccentricity vector from the velocity and position vectors
/// 
/// Inputs
/// ------
/// pos: `Vector3<f64>`           
///     Position vector        
///    
/// vel: `Vector3<f64>`
///     Velocity vector     
///       
/// grav_param: `f64`
///     Gravitational parameter of the central body
/// 
/// Outputs
/// -------
/// ecc_vec: `Vector3<f64>`           
///     Eccentricity vector
fn calc_ecc_vec(
    pos: Vector3<f64>,
    vel: Vector3<f64>, 
    grav_param: f64
) -> Vector3<f64> {
    let spec_linear_moment: f64 = pos.dot(&vel);
    let v_sq: f64 = vel.norm().powi(2);
    let ecc_vec: Vector3<f64> = ((
        v_sq - grav_param / pos.norm())*pos - (spec_linear_moment*vel)) / grav_param;
    ecc_vec
}

/// Calculate the RAAN given the ascending node vector.
/// 
/// Inputs
/// ------
/// ascend_node_vec: `Vector3<f64>`
///     Vector defining the ascending node.
/// 
/// Outputs
/// -------
/// raan: `f64`
///     Right ascension of the ascending node.
pub fn calc_raan(ascend_node_vec: Vector3<f64>) -> f64 {
    let raan: f64 = (ascend_node_vec[0] / ascend_node_vec.norm()).acos();
    raan
}

/// Calculate inclination
/// 
/// # Arguments
/// 
/// * `spec_ang_moment` - specific angular momentum vector
/// 
/// # Returns
/// 
/// * `f64` - inclination in radians
pub fn calc_inclination(spec_ang_moment: Vector3<f64>) -> f64 {
    let inclination: f64 = (spec_ang_moment[2] / spec_ang_moment.norm()).acos();
    inclination
}

/// Calculate the mean motion of an orbit.
/// 
/// # Arguments
/// 
/// * `semi_major_axis` - The semi-major axis of the orbit.
/// * `grav_param` - The gravitational parameter of the body.
pub fn calc_mean_motion(semi_major_axis: f64, grav_param: f64) -> f64 {
    let mean_motion: f64 = 1.0 / (2.0 * PI * (semi_major_axis.powi(3)/grav_param).sqrt());
    mean_motion
}

/// Calculate semi major axis of orbit
/// 
/// Inputs
/// ------
pub fn calc_semi_major_axis(grav_param: f64, mean_motion: f64) -> f64 {
    let semi_major_axis: f64 = ((grav_param)/mean_motion.powi(2)).powf(1.0/3.0);
    semi_major_axis
}

/// Calculate total delta v for hohmann transfer
/// 
/// Inputs
/// ------
/// radius_1 : `f64`
///     Radius of inner orbit
/// 
/// radius_2 : `f64`
///     Radius of outer orbit
/// 
/// vel_0 : `f64`
///     Initial Velocity of vehicle
/// 
/// Outputs
/// -------
/// delta_v_total : `f64`
///     Total delta v for maneuver
pub fn calc_hohmann_transfer(
    radius_1: f64,
    radius_2: f64,
    vel_0: f64
) -> f64 {
    let delta_v_1: f64 = vel_0 * (((2.0 * radius_2)/(radius_1 + radius_2)).sqrt()- 1.0);
    let delta_v_2_num: f64 = (radius_1 / radius_2).sqrt() + (2.0 * radius_1);
    let delta_v_2_den: f64 = (radius_2 *(1.0 + radius_2 / radius_1)).sqrt();
    let delta_v_2: f64 = vel_0 * delta_v_2_num / delta_v_2_den; 
    let delta_v_total: f64 = delta_v_1 + delta_v_2;
    return delta_v_total
}

/// Calculate sphere of influence for a body
/// 
/// Inputs
/// ------
/// mass_0 : `f64`
///     Primary mass of 2 body system
/// 
/// mass_1 : `f64`
///     Secondary mass of 2 body system
/// 
/// semi_major_axis : `f64`
///     Orbital ellipse semi major axis
/// 
/// eccentricity : `f64`
///     Orbital ellipse eccentricity
/// 
/// Outputs
/// -------
/// radius : `f64`
///     Radius of sphere of influence
pub fn calc_hill_sphere(
    mass_0: f64,
    mass_1: f64,
    semi_major_axis: f64,
    eccentricity: f64
) -> f64 {
    let radius: f64 = 
        semi_major_axis * (1.0 - eccentricity) * (mass_0 / (3.0 * mass_1)).powf(1.0/3.0);
    return radius
}


#[cfg(test)]
mod orbit_tests {
    use super::*;

    #[test]
    fn test_hohmann_transfer(){
        let radius_1: f64 = 5000.;
        let radius_2: f64 = 5500.;
        let vel_0: f64 = 10.;
        let deltav: f64 = calc_hohmann_transfer(radius_1, radius_2, vel_0);

        assert_eq!(deltav, 930.8082549013038);

    }

    #[test]
    fn test_hill_sphere(){
        let earth_mass: f64 = crate::constants::EARTH_MASS;
        let sun_mass: f64 = crate::constants::SUN_MASS;
        let earth_orbit_semi_major: f64 = crate::constants::EARTH_ORBIT_SEMI_MAJOR;
        let earth_orbit_ecc: f64 = crate::constants::EARTH_ORBIT_ECC;

        let sphere_rad: f64 = calc_hill_sphere(
            sun_mass, 
            earth_mass,
            earth_orbit_semi_major,
            earth_orbit_ecc
        );

        assert_eq!(sphere_rad, 7069993624.788241)
    }


    #[test]
    fn test_kepler_orbit(){
        let tle_str = "ISS (ZARYA) 
        1 25544U 98067A   23035.69666365  .00008902  00000+0  16600-3 0  9994
        2 25544  51.6420 264.7747 0008620 314.4274 150.8239 15.49588766381243";
        let kep = Orbit::from_tle(tle_str.to_string());

        assert_eq!(kep.semi_major_axis, 11840.341648011852);
        assert_eq!(kep.eccentricity, 0.0008620);
        assert_eq!(kep.inclination, 51.6420);
        assert_eq!(kep.raan, 264.7747);
        assert_eq!(kep.argument_of_perigee, 314.4274);
    }
}
