/*
Gravitational Bodies
*/

use nalgebra as na;
use na::Vector3;
use chrono::{DateTime as DateTime, NaiveDateTime, NaiveDate, NaiveTime, Utc};
use std::f64::consts::PI;

use crate::constants as cst;
use crate::cosmos::time as time;

/// Gravitational Body
#[derive(Clone, Debug, PartialEq)]
pub struct Body{
    pub name: String,
    pub grav_param: f64,
    pub eq_radius: f64,
    pub rotation_rate: f64,
    pub eccentricity: f64
}


/// Gravitational Body
impl Body {

    /// Convert into simple particle type
    /// 
    /// Inputs
    /// ------
    /// motion: `Vec<Vector3<f64>>`
    ///     motion of body in frame
    /// 
    /// Outputs
    /// -------
    /// particle: `Particle`
    ///     `Particle` equivalent of body
    pub fn to_particle(&self, motion: Vec<Vector3<f64>>) -> Particle {
        let mass: f64 = self.grav_param / cst::GRAV_CONST;
        let particle: Particle = Particle {mass: mass, motion: motion};
        particle
    }


    /// Calculate gravitational acceleration at radial distance from Body
    /// 
    /// Inputs
    /// ------
    /// radius: `Vector3<f64>`
    ///     Radius in km from Body center
    /// 
    /// Outputs
    /// -------
    /// grav_acc: `Vector3<f64>`
    ///     Acceleration rate due to gravity
    pub fn calc_body_grav(&self, radius: Vector3<f64>) -> Vector3<f64> {
        let grav_acc: Vector3<f64>  = self.grav_param * radius / radius.norm().powi(3);
        return grav_acc
    }
    

    /// Calculate required orbital velocity at radial distance
    /// 
    /// Inputs
    /// ------
    /// radius: `f64`
    ///     Radius in km from Body center
    ///     
    /// Outputs
    /// -------
    /// val: `f64`
    ///     Required tangential velocity magnitude
    pub fn calc_orbital_velocity_mag(&self, radius: f64) -> f64 {
        let vel: f64 = (2.0 * self.grav_param / radius).sqrt();
        return vel
    }

    /// Calculate period of orbit at semi major axis
    /// 
    /// Inputs
    /// ------
    /// semi_major_axis: `f64`
    ///     Semi major axis of orbital ellipse
    /// 
    /// Outputs
    /// -------
    /// period: `f64`
    ///    Period of orbit in seconds
    pub fn calc_period(&self, semi_major_axis: f64) -> f64 {
        let period: f64 = 2.0 * PI * (semi_major_axis.powi(3)/self.grav_param).sqrt();
        return period
    }
    
    /// Calculate radius for stationary orbit above body surface
    /// 
    /// Outputs
    /// -------
    /// r_mag: `f64`
    ///     Magnitude of radius for stationary orbit
    pub fn calc_stationary_orbit(&self) -> f64 {
        let period: f64 = 2.0 * PI / self.rotation_rate;
        let a: f64 = self.grav_param * period.powi(2); 
        let r_mag: f64 = (a / (4.0 * PI.powi(2))).powf(1.0 / 3.0); 
        return r_mag
    }

    /// Geodetic to rectangular coordinates
    /// E.g. Latitude, Longitude, Altitude to ECEF
    /// 
    /// Inputs
    /// ------
    /// lla: `Vector3<f64>`
    ///     geodetic coords
    /// 
    /// Outputs
    /// -------
    /// xyz: `Vector3<f64>`
    ///     Cartesian coords
    pub fn geodetic_to_xyz(&self, lla: Vector3<f64>) -> Vector3<f64> {
        let radius: f64 = self.calc_prime_vertical(lla[0]);
        let x: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].cos();
        let y: f64 = (radius + lla[2]) * lla[0].cos() * lla[1].sin();
        let z: f64 = ((1.0 - self.eccentricity.powi(2)) * radius + lla[2]) * lla[0].sin();
        let xyz: Vector3<f64> = Vector3::new(x, y, z); 
        return xyz
    }

    /// Calculate prime vertical radius to surface at latitude
    /// 
    /// Inputs
    /// ------
    /// lat_deg: `f64`
    ///     Lattitude in degrees
    pub fn calc_prime_vertical(&self, lat_deg: f64) -> f64 {
        let lat_radians: f64 = PI * lat_deg / 180.0;
        let radius: f64 = 
            self.eq_radius / (1.0 - (self.eccentricity * lat_radians.sin()).powi(2)).sqrt();
        return radius
    }

    /// Rectangular coordinates to geodetic
    /// E.g. ECEF to LLH
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
        let cap_p: f64 = 54.0 * b.powi(2)*xyz[2].powi(2) / 
            (3.0 * (s + 1.0 + 1.0 / s).powi(2) * g.powi(2));
        let q: f64 = (1.0 + 2.0 * ecc_2.powi(2) * cap_p).sqrt();
        let r_0: f64 = -cap_p * ecc_2 * p /(1.0+q) + 
            ((a.powi(2)/2.0)*(1.0 + 1.0 / q) - 
            cap_p * (1.0 - ecc_2) * xyz[2].powi(2) / (q * (1.0 + q)) - 
            cap_p*p.powi(2)/2.0).sqrt();
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


/// Gravitational Particle
#[derive(Clone, Debug, PartialEq)]
pub struct Particle {
    pub mass: f64,
    pub motion: Vec<Vector3<f64>>,
}

impl Particle {

    /// Convert into Body type
    /// 
    /// Inputs
    /// ------
    /// name : `String`
    ///     Name of body
    /// 
    /// eq_radius : `f64`
    ///     Equatorial radius of body
    /// 
    /// rotation_rate : `f64`
    ///     Rotation rate of body
    /// 
    /// eccentricity : `f64`
    ///     Body oblateness
    pub fn to_body(&self, name: String, eq_radius: f64, rotation_rate: f64, eccentricity: f64) -> Body {
        let grav_param: f64 = self.mass * cst::GRAV_CONST;
        let body: Body = Body {
            name,
            grav_param,
            eq_radius,
            rotation_rate,
            eccentricity
        };
        body
    }
}


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
    /// name : `String`
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

        let md: (u32, u32) = time::calc_month_day(day_of_year, year);
        let date: NaiveDate = NaiveDate::from_ymd_opt(
            year, md.0, md.1).unwrap();
        
        let percent_of_day: f64 = epoch_day_full[0]
            .to_string()
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

        let ecc: f64 = line2[4]
            .to_string()
            .parse::<f64>()
            .unwrap() * 10e-7;

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
    /// name : `String`
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
    /// 
    /// # Example
    /// ```rust,no_run      
    /// use gravity::kepler::Orbit;
    /// 
    /// let name: &str = "Moon"
    /// let grav_param: f64 = cst::MOON_GRAV_PARAM
    /// let pos: Vector3<f64> = Vector3::new(0.0, 0
    ///     0.0, 0.0)
    /// let vel: Vector3<f64> = Vector3::new(0.0, 0
    ///     0.0, 0.0) }
    /// let epoch_datetime: DateTime<Utc> = Utc.ymd(1969
    ///     7, 20).and_hms(0, 0, 0) }
    /// let orbit = Orbit::from_pos_vel(name, grav_param, pos, vel
    ///     epoch_datetime)
    /// 
    /// assert_eq!(orbit.name, "Moon")
    /// assert_eq!(orbit.grav_param, cst::MOON_GRAV_PARAM)
    /// assert_eq!(orbit.semi_major_axis, cst::MOON_SEMI_MAJOR_AXIS)
    pub fn from_pos_vel(
        name: String,
        grav_param: f64,
        pos: Vector3<f64>,
        vel: Vector3<f64>,
        epoch_datetime: DateTime<Utc>
    ) -> Self {
        let spec_ang_moment: Vector3<f64> = pos.cross(&vel);
        let spec_lin_moment: f64 = pos.dot(&vel);

        let ecc_vec: Vector3<f64> = 
            ((vel.norm().powi(2) - grav_param / pos.norm())*pos 
            - (spec_lin_moment*vel)) / grav_param;
        let ascend_node_vec: Vector3<f64> = Vector3::z_axis().cross(&spec_ang_moment);

        let semi_major_axis: f64 = 
            spec_ang_moment.norm().powi(2) * 
            (1.0 - ecc_vec.norm_squared()) / grav_param;

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



    pub fn propogate(
        &self, 
        eval_datetimes: Vec<DateTime<Utc>, 
        precision: f64)
    -> Vec<Vector3<f64>> {
        let arg_perigree = self.argument_of_perigee;
        return 
    }

}


pub fn calc_raan(ascend_node_vec: Vector3<f64>) -> f64 {
    (ascend_node_vec[0] / ascend_node_vec.norm()).acos()
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
    (spec_ang_moment[2] / spec_ang_moment.norm()).acos()
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
    fn test_orbit(){
        let tle_str = "ISS (ZARYA) 
        1 25544U 98067A   23035.69666365  .00008902  00000+0  16600-3 0  9994
        2 25544  51.6420 264.7747 0008620 314.4274 150.8239 15.49588766381243";
        let orb = Orbit::from_tle(tle_str.to_string());

        assert_eq!(orb.argument_of_perigee, 314.4274);
    }


}
