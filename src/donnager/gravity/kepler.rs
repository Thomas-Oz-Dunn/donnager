/*
Gravitational Bodies
kraus
*/

use nalgebra::{Vector3, Matrix3};
use chrono::{DateTime, TimeZone, Utc};
use plotters::prelude::*;
use std::{f64::consts::PI, vec, ops::Range};

use crate::donnager::{spacetime as xyzt, constants as cst};

pub struct Maneuver{
    pub delta_v: Vector3<f64>,
    pub act_time: f64
}

/// Calculate coplanar maneuver
/// 
/// Inputs
/// ------
/// orbit_0
/// 
/// orbit_f
/// 
/// 
pub fn calc_coplanar_maneuver(
    orbit_0: Orbit,
    orbit_f: Orbit,
    epoch_datetime: DateTime<Utc>
) -> Maneuver {

    // let cos_i_0 = (orbit_0.inclination * cst::DEG_TO_RAD).cos();
    // let cos_i_1 = (orbit_f.inclination * cst::DEG_TO_RAD).cos();
    // let sin_i_0 = (orbit_0.inclination * cst::DEG_TO_RAD).sin();
    // let sin_i_1 = (orbit_f.inclination * cst::DEG_TO_RAD).sin();
    // let theta= ((
    //     orbit_f.raan * cst::DEG_TO_RAD - orbit_0.raan * cst::DEG_TO_RAD).cos()* sin_i_1 * sin_i_0 + cos_i_0*cos_i_1).acos();

    let frame: xyzt::ReferenceFrames = xyzt::ReferenceFrames::InertialCartesian;
    let t_start: f64 = epoch_datetime.timestamp() as f64;

    let motion0 = orbit_0.calc_motion(t_start, frame);
    let motionf = orbit_f.calc_motion(t_start, frame);

    let delv = motionf[1] - motion0[1];
    let man: Maneuver = Maneuver { delta_v: (delv), act_time: (t_start) };

    return man
}


/// Calculate maneuvers
/// 
/// Inputs
/// ------
/// orbit_0
/// 
/// orbit_f
/// 
/// epoch_datetime
pub fn calc_maneuvers(
    orbit_0: Orbit,
    orbit_f: Orbit,
    epoch_datetime: DateTime<Utc>
) -> Vec<Maneuver> {
    let frame: xyzt::ReferenceFrames = xyzt::ReferenceFrames::InertialCartesian;

    if orbit_0.inclination != orbit_f.inclination {
        // Make coplanar
        let delv = calc_coplanar_maneuver(orbit_0, orbit_f, epoch_datetime);
        return vec![delv]

    }

    // Hohmann transfer
    let radius_1: f64 = orbit_0.semi_major_axis;
    let radius_2: f64 = orbit_f.semi_major_axis;
    let t_start: f64 = epoch_datetime.timestamp() as f64;
    let motion1 = orbit_0.calc_motion(t_start, frame);
    let (del1, del2) = calc_hohmann_transfer(radius_1, radius_2, motion1[1].norm());

    // Maneuver 1
    let dv_1 = del1 * motion1[1] / motion1[1].norm();
    let man1: Maneuver = Maneuver { delta_v: (dv_1), act_time: (t_start) };
    
    let name = "Transfer orbit";
    let vel = motion1[1] + dv_1;
    let trans_orbit: Orbit = Orbit::from_pos_vel(
        name.to_string(), 
        orbit_0.central_body, 
        motion1[0], 
        vel, 
        epoch_datetime);

    let period: f64 = trans_orbit.calc_period();

    // Maneuver 2
    let time_2: f64 = t_start + period / 2.;
    let motion2 = trans_orbit.calc_motion(time_2, frame);
    let dv_2 = del2 * motion2[1] / motion2[1].norm();
    let man2: Maneuver = Maneuver { delta_v: (dv_2), act_time: (time_2) };

    let maneuvers = vec![man1, man2];

    return maneuvers

}


/// Orbit structure
/// 
/// Variables
/// ---------
/// name : `String`
///     Name of orbiting body
/// 
/// central_body : `xyzt::Body`
///     Central gravitational body
/// 
/// semi_major_axis : `f64`
/// 
/// eccentricity : `f64`
/// 
/// inclination : `f64`
/// 
/// raan : `f64`
/// 
/// argument_of_perigee : `f64`
/// 
/// mean_anomaly : `f64`
/// 
/// mean_motion : `f64`
/// 
/// epoch : `DateTime<Utc>`
///     Reference epoch DateTime
/// 
/// Methods
/// -------
/// from_keplerian : `fn`
///     Populate Orbit from Keplerian parameters
/// 
/// propagate : `fn`
///     Propagate Orbit forward in time
#[derive(Clone, Debug, PartialEq)]
pub struct Orbit{
    pub name: String,
    pub central_body: xyzt::Body,
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
    ///     Name of orbit
    /// 
    /// central_body : `xyzt::Body`           
    ///     Central body of orbit
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
        central_body: xyzt::Body,
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
            name: name.to_string(),
            central_body,
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

        let md: (u32, u32) = xyzt::calc_month_day(day_of_year, year);
        
        let percent_of_day: f64 = 
        (".".to_owned() + epoch_day_full[1])
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

        let epoch_date_time = xyzt::ymd_hms_to_datetime(
            year, md.0, md.1, hours_whole, minutes_whole, seconds_whole);
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
            (".".to_owned() + line2[4])
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
            cst::EARTH::GRAV_PARAM, mean_motion);

        // Earth
        let earth: xyzt::Body = xyzt::Body {
            name: String::from("Earth"),
            grav_param: cst::EARTH::GRAV_PARAM,
            eq_radius: cst::EARTH::RADIUS_EQUATOR,
            rotation_rate: cst::EARTH::ROT_RATE,
            sidereal_day_hours: cst::EARTH::SIDEREAL_DAY,
            eccentricity: cst::EARTH::SURFACE_ECC
        };

        Orbit {
            name: name.to_string(),
            central_body: earth,
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
        central_body: xyzt::Body,
        pos: Vector3<f64>,
        vel: Vector3<f64>,
        epoch_datetime: DateTime<Utc>
    ) -> Self {
        let grav_param: f64 = central_body.grav_param;
        let spec_ang_moment: Vector3<f64> = pos.cross(&vel);
        let ecc_vec: Vector3<f64> = calc_ecc_vec(pos, vel, grav_param);
        let ascend_node_vec: Vector3<f64> = Vector3::z_axis().cross(&spec_ang_moment);

        let semi_major_axis: f64 = 
            spec_ang_moment.norm_squared() * (1.0 - ecc_vec.norm_squared()) / grav_param;

        Orbit {
            name: name.to_string(),
            central_body,
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
    /// 
    /// Outputs
    /// -------
    /// motion: `[order, xyz]`
    ///     Position and Velocity in reference frame
    pub fn calc_motion(
        &self, 
        time: f64, 
        frame: xyzt::ReferenceFrames
    ) -> Vec<Vector3<f64>> {
        let true_anomaly_rad: f64 = self.calc_true_anomaly(time);
        let cos_true_anomaly: f64 = true_anomaly_rad.cos();
        let sin_true_anomaly: f64 = true_anomaly_rad.sin();
        let radius: f64 = self.calc_radius(cos_true_anomaly);
        
        // Vectorize with [n_evals, 3dim, n_order]
        // Perifocal
        let x_pos: f64 = radius * cos_true_anomaly;
        let y_pos: f64 = radius * sin_true_anomaly;
        let z_pos: f64 = 0.0;
        let pos: Vector3<f64> = Vector3::new(x_pos, y_pos, z_pos);

        // Perifocal
        let x_vel: f64 = -self.mean_motion * radius * sin_true_anomaly;
        let y_vel: f64 = self.mean_motion * radius * (self.eccentricity + cos_true_anomaly);
        let z_vel: f64 = 0.0;
        let vel: Vector3<f64> = Vector3::new(x_vel, y_vel, z_vel);
        
        match frame {
            xyzt::ReferenceFrames::Perifocal => {
                return vec![pos, vel];
            },
            xyzt::ReferenceFrames::InertialCartesian => {
                let rotam: Matrix3<f64> = self.calc_pfcl_inertial_rotam();
                let p_ic: Vector3<f64> = rotam * pos;
                let v_ic: Vector3<f64> = rotam * vel;
                return vec![p_ic, v_ic];
            },
            xyzt::ReferenceFrames::RotationalCartesian => {
                let pfcl_eci_rotam: Matrix3<f64> = self.calc_pfcl_inertial_rotam();
                let eci_pos: Vector3<f64> = pfcl_eci_rotam * pos;
                let eci_vel: Vector3<f64> = pfcl_eci_rotam * vel;

                let new_time: f64 = self.epoch.timestamp() as f64 + time;
                let new_epoch_datetime: DateTime<Utc> = Utc.timestamp_opt(
                    new_time as i64, 0).unwrap();
                let eci_ecef_rotam: Matrix3<f64> = 
                    xyzt::calc_inertial_rotational_rotam(
                        new_epoch_datetime, 
                        self.central_body.clone());

                let ecef_pos: Vector3<f64> = eci_ecef_rotam * eci_pos;
                let ecef_vel: Vector3<f64> = eci_ecef_rotam * eci_vel;
                return vec![ecef_pos, ecef_vel];
            },
            xyzt::ReferenceFrames::Planetodetic => {
                let pfcl_eci_rotam: Matrix3<f64> = self.calc_pfcl_inertial_rotam();
                let eci_pos: Vector3<f64> = pfcl_eci_rotam * pos;
                let eci_vel: Vector3<f64> = pfcl_eci_rotam * vel;

                let new_time: f64 = self.epoch.timestamp() as f64 + time;
                let new_epoch_datetime: DateTime<Utc> = Utc.timestamp_opt(
                    new_time as i64, 0).unwrap();
                let eci_ecef_rotam: Matrix3<f64> = 
                    xyzt::calc_inertial_rotational_rotam(
                        new_epoch_datetime, 
                        self.central_body.clone());

                let ecef_pos: Vector3<f64> = eci_ecef_rotam * eci_pos;
                let ecef_vel: Vector3<f64> = eci_ecef_rotam * eci_vel;

                let lla_pos: Vector3<f64> = xyzt::ecef_to_lla(ecef_pos);
                let lla_vel: Vector3<f64> = xyzt::ecef_to_lla(ecef_vel);
                return vec![lla_pos, lla_vel];
            }
        }
    }

    /// Calculate radius at true anomaly
    /// 
    /// Inputs
    /// ------
    /// cos_true_anomaly: `f64`
    ///     Cosine of true anomaly angle
    /// 
    /// Outputs
    /// -------
    /// radius: `f64`
    ///     Magnitude of radius
    fn calc_radius(&self, cos_true_anomaly: f64) -> f64 {
        let radius: f64 = self.semi_major_axis * 
            (1.0 - self.eccentricity.powi(2)) / 
            (1.0 + self.eccentricity * cos_true_anomaly);
        return radius
    }

    /// Calculate true anomaly in radians
    /// 
    /// Inputs
    /// ------
    /// time: `f64`
    ///     Time since epoch
    pub fn calc_true_anomaly(&self, time: f64) -> f64 {
        let mean_anom: f64 = (
            self.mean_anomaly + self.mean_motion * time) * cst::DEG_TO_RAD;
        let ecc_anom: f64 = (
            mean_anom - self.eccentricity * cst::DEG_TO_RAD * (
                1.0 - mean_anom.cos())) * cst::DEG_TO_RAD;
        
        let true_anomaly_rad: f64 = 2.0 * (
            ecc_anom.sin()).atan2(-ecc_anom.cos());
        return true_anomaly_rad
    }

    /// Calculate perifocal to eci rotation matrix
    pub fn calc_pfcl_inertial_rotam(&self) -> Matrix3<f64> {
        let cos_raan: f64 = self.raan.cos();
        let sin_raan: f64 = self.raan.sin();
        let cos_inc: f64 = self.inclination.cos();
        let sin_inc: f64 = self.inclination.sin();
        let cos_arg_peri: f64 = self.argument_of_perigee.cos();
        let sin_arg_peri: f64 = self.argument_of_perigee.sin();

        let rot_mat: Matrix3<f64> = 
            Matrix3::new(
                cos_raan, -sin_raan, 0.0,
                sin_raan, cos_raan, 0.0,
                0.0, 0.0, 1.0);

        let rot_mat_2: Matrix3<f64> =
            Matrix3::new(
                cos_inc, 0.0, sin_inc,
                0.0, 1.0, 0.0, 
                -sin_inc, 0.0, cos_inc);

        let rot_mat_3: Matrix3<f64> = 
            Matrix3::new(
                cos_arg_peri, sin_arg_peri, 0.0, 
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
        let frame = xyzt::ReferenceFrames::InertialCartesian;
        let motion = self.calc_motion(new_time, frame);
        let new_epoch_datetime: DateTime<Utc> = Utc.timestamp_opt(
            new_time as i64, 0).unwrap();

        let new_orbit: Orbit = Orbit::from_pos_vel(
            self.name.clone(), 
            self.central_body.clone(), 
            motion[0], 
            motion[1], 
            new_epoch_datetime);
        *self = new_orbit;

        }

    /// Show orbit plot
    pub fn show(&self, frame: xyzt::ReferenceFrames) {

        let pathname = format!("{}_orbit_{:?}.png", self.name, &frame);
        let plottitle = format!("{} orbit {:?} frame", self.name, &frame);
        
        let t_start = self.epoch.timestamp() as f64;
        let t_step = 1.;
        let times = (
            (t_start)..(t_start + self.mean_motion * 43.)
        ).step(t_step);
        
        let mut pos_mut: Vec<Vector3<f64>> = Vec::new();
        times.values().for_each(|time| {
            let motion_frame = self.calc_motion(time, frame);
            pos_mut.push(motion_frame[0]);
        });

        // Plot
        let drawing_area = 
            BitMapBackend::new(&pathname, (500, 300))
            .into_drawing_area();

        let mut chart_builder = ChartBuilder::on(&drawing_area);

        drawing_area
            .fill(&WHITE)
            .unwrap();

        chart_builder
            .margin(5)
            .set_left_and_bottom_label_area_size(35)
            .caption(
                plottitle, (
                    "Times New Roman", 
                    20, 
                    FontStyle::Bold, 
                    &BLACK)
                    .into_text_style(&drawing_area));

        let body_rad = self.central_body.eq_radius;

        match frame {
            // xyzt::ReferenceFrames::Heliocentric => {
                
            //     let x_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad; 
            //     let y_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad;
            //     let z_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad;

            //     let mut chart = chart_builder
            //         .build_cartesian_3d(x_spec, y_spec, z_spec)
            //         .unwrap();

            //     chart
            //         .configure_axes()
            //         .draw()
            //         .unwrap();

            //     chart
            //         .configure_series_labels()
            //         .background_style(&WHITE.mix(0.8))
            //         .border_style(&BLACK)
            //         .draw()
            //         .unwrap();

            //     chart
            //         .draw_series(
            //             SurfaceSeries::xoy(
            //                 (-100..100).map(|f| {(f as f64 * body_rad / 100.).cos()}),
            //                 (-100..100).map(|f| {(f as f64 * body_rad / 100.).sin()}),
            //                 |x, y| (-(x * x + y * y).sqrt()),
            //             )
            //         ).unwrap();

            //     chart
            //         .draw_series(
            //             SurfaceSeries::xoy(
            //                 (-100..100).map(|f| {(f as f64 * body_rad / 100.).cos()}),
            //                 (-100..100).map(|f| {(f as f64 * body_rad / 100.).sin()}),
            //                 |x, y| (-(x * x + y * y).sqrt()),
            //             )).unwrap()
            //         .label("Sun");

            // },
            xyzt::ReferenceFrames::InertialCartesian => {
                
                let x_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad; 
                let y_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad;
                let z_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad;

                let mut chart = chart_builder
                    .build_cartesian_3d(x_spec, y_spec, z_spec)
                    .unwrap();

                chart
                    .configure_axes()
                    .draw()
                    .unwrap();

                chart
                    .configure_series_labels()
                    .background_style(&WHITE.mix(0.8))
                    .border_style(&BLACK)
                    .draw()
                    .unwrap();

                chart
                    .draw_series(
                        SurfaceSeries::xoy(
                            (-100..100).map(|f| {(f as f64 * body_rad / 100.).cos()}),
                            (-100..100).map(|f| {(f as f64 * body_rad / 100.).sin()}),
                            |x, y| (-(x * x + y * y).sqrt()),
                        )
                    ).unwrap();

                chart
                    .draw_series(
                        SurfaceSeries::xoy(
                            (-100..100).map(|f| {(f as f64 * body_rad / 100.).cos()}),
                            (-100..100).map(|f| {(f as f64 * body_rad / 100.).sin()}),
                            |x, y| (-(x * x + y * y).sqrt()),
                        )).unwrap()
                    .label("Earth");

            },
            xyzt::ReferenceFrames::RotationalCartesian => {
                // 3D globe w/ spin
                let x_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad; 
                let y_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad;
                let z_spec: Range<f64> = -1.5 * body_rad..1.5 * body_rad;

                let mut chart = 
                    chart_builder.build_cartesian_3d(x_spec, y_spec, z_spec).unwrap();

                chart.configure_axes().draw().unwrap();

                chart
                    .configure_series_labels()
                    .background_style(&WHITE.mix(0.8))
                    .border_style(&BLACK)
                    .draw().unwrap();


                chart
                    .draw_series(
                        SurfaceSeries::xoy(
                            (-100..100).map(|f| {(f as f64 * body_rad / 100.).cos()}),
                            (-100..100).map(|f| {(f as f64 * body_rad / 100.).sin()}),
                            |x, z| (-(x * x + z * z).sqrt()))
                            .style(&BLUE.mix(0.5))  
                    ).unwrap();

                chart
                    .draw_series(
                        SurfaceSeries::xoy(
                            (-100..100).map(|f| {(f as f64 * body_rad / 100.).cos()}),
                            (-100..100).map(|f| {(f as f64 * body_rad / 100.).sin()}),
                            |x, z| (-(x * x + z * z).sqrt()))
                        .style(&BLUE.mix(0.5))  
                        ).unwrap()
                    .label("Earth");
    
            },
            xyzt::ReferenceFrames::Planetodetic => {
                // 2d Ground track
                // TODO-TD: add global shoreline trace
                let y_spec: Range<f64> = -90.0..90.;  // N to S poles
                let x_spec: Range<f64> = -180.0..180.;  // International date line
                let mut chart = 
                    chart_builder.build_cartesian_2d(x_spec, y_spec).unwrap();

                chart.configure_mesh().draw().unwrap();

                chart.draw_series(
                    PointSeries::of_element(
                        pos_mut.iter().map(|p| (p.y, p.x)),
                        1,
                        &BLUE,
                        &|c, s, st| {
                            Circle::new((c.0, c.1), s, st.filled())}
                    )
                ).unwrap()
                .label("Orbit")
                .legend(
                    |(x, y)| 
                    PathElement::new(vec![(x, y), (x + 20, y)], 
                    &BLUE));


                chart.configure_series_labels()
                    .background_style(&WHITE.mix(0.8))
                    .border_style(&BLACK)
                    .draw().unwrap();

            },
            xyzt::ReferenceFrames::Perifocal => {
                // 2d planar plot
                let semi_major: f64 = self.semi_major_axis;
                let semi_latus: f64 = semi_major * (1. - self.eccentricity);
                let x_spec: Range<f64> = -semi_major..semi_major; 
                let y_spec: Range<f64> = -semi_latus..semi_latus;
                let mut chart = 
                    chart_builder.build_cartesian_2d(x_spec, y_spec).unwrap();

                chart.configure_mesh().draw().unwrap();

                chart.configure_series_labels()
                    .background_style(&WHITE.mix(0.8))
                    .border_style(&BLACK)
                    .draw().unwrap();
                
                chart.draw_series(
                    PointSeries::of_element(
                        pos_mut.iter().map(|p| (p.x, p.y)),
                        1,
                        &BLACK,
                        &|c, s, st| {
                            Circle::new((c.0, c.1), s, st.filled())}
                    )
                ).unwrap();

                let mut bod_surf: Vec<(f64, f64)> = Vec::new();
                for theta in 0..360{
                    bod_surf.push(
                        (body_rad * (theta as f64 * cst::DEG_TO_RAD).cos(),
                        body_rad * (theta as f64 * cst::DEG_TO_RAD).sin()))
                }

                chart.draw_series(
                    PointSeries::of_element(
                        bod_surf.iter().map(|p|{(p.0, p.1)}), 
                        1,
                        &BLUE,
                        &|c, s, st| {
                            Circle::new((c.0, c.1), s, st.filled())}
                    )
                ).unwrap();

            }   
        }
    }


    /// Convert orbit to KML string
    /// 
    /// Inputs
    /// ------
    /// None
    /// 
    /// Outputs
    /// -------
    /// kml_string: `String`
    ///     KML string of orbit
    pub fn to_kml(&self) -> String {
        let mut kml_string_mut: String = String::new();
        kml_string_mut.push_str("<Placemark>\n");
        kml_string_mut.push_str("<name>");
        kml_string_mut.push_str(&self.name);
        kml_string_mut.push_str("</name>\n");
        kml_string_mut.push_str("<LineString>\n");
        kml_string_mut.push_str("<extrude>1</extrude>\n");
        kml_string_mut.push_str("<tessellate>1</tessellate>\n");
        kml_string_mut.push_str("<altitudeMode>absolute</altitudeMode>\n");
        kml_string_mut.push_str("<coordinates>\n");
        let mut time: f64 = self.epoch.timestamp() as f64;
        let frame = xyzt::ReferenceFrames::RotationalCartesian;

        let mut motion_ecef = self.calc_motion(time, frame);
        let pos_lla = xyzt::ecef_to_lla(motion_ecef[0]);
        kml_string_mut.push_str(&format!("{},{},{}\n", pos_lla[1], pos_lla[0], pos_lla[2]));
        time += 60.0;
        while time < self.epoch.timestamp() as f64 + 86400.0 {
            motion_ecef = self.calc_motion(time, frame);
            let pos_lla = xyzt::ecef_to_lla(motion_ecef[0]);
            kml_string_mut.push_str(&format!("{},{},{}\n", pos_lla[1], pos_lla[0], pos_lla[2]));
            time += 60.0;
        }
        kml_string_mut.push_str("</coordinates>\n");
        kml_string_mut.push_str("</LineString>\n");
        kml_string_mut.push_str("</Placemark>\n");
        let kml_string = kml_string_mut;
        return kml_string
    }


    /// Calculate ground coverage
    /// 
    /// Inputs
    /// ------
    /// time: `f64`
    ///     Time for evaluation
    pub fn calc_ground_coverage_radius(&self, time: f64) -> f64 {
        let motion_ecef = self.calc_motion(
            time, xyzt::ReferenceFrames::RotationalCartesian);
        let pos_lla: Vector3<f64> = xyzt::ecef_to_lla(motion_ecef[0]);
        let theta: f64 = (
            self.central_body.eq_radius / (self.central_body.eq_radius + pos_lla.z)).acos();
        let cov_radius: f64 = theta * self.central_body.eq_radius;
        return cov_radius
    }


    /// Calculate orbital period
    pub fn calc_period(&self) -> f64 {
        return calc_period(self.semi_major_axis, self.central_body.grav_param)
    }

}


/// Calculate lagrange point locations in DU
/// 
/// Inputs
/// ------
/// mass_1: `f64`
///     Larger mass
/// 
/// mass_2: `f64`
///     Smaller mass
pub fn calc_lagrange_points(
    mass_1: f64,
    mass_2: f64
) -> Vec<Vector3<f64>> {
    let mass_ratio: f64= mass_2 / (mass_1 + mass_2);

    // Permissible if mass_2 is <<<< mass_1
    let xl12: f64 = (mass_2/(3.*mass_1)).powf(1./3.);
    let xl3: f64 = 7.*mass_2 / (12. * mass_1);

    let l1: Vector3<f64> = Vector3::new(1. - xl12, 0.,0.);
    let l2: Vector3<f64> = Vector3::new(1. + xl12 ,0.,0.);
    let l3: Vector3<f64> = Vector3::new(-xl3 ,0.,0.);
    let l4: Vector3<f64> = Vector3::new(mass_ratio - 0.5, -(3_f64).sqrt()/2., 0.);
    let l5: Vector3<f64> = Vector3::new(mass_ratio - 0.5, (3_f64).sqrt()/2., 0.);
    return vec![l1, l2, l3, l4, l5]
}


/// Calculate period of orbit at semi major axis
/// 
/// Inputs
/// ------
/// semi_major_axis: `f64`
///     Semi major axis of orbital ellipse
/// 
/// grav_param: `f64`
///     Gravitational Parameter
/// 
/// Outputs
/// -------
/// period: `f64`
///    Period of orbit in seconds
pub fn calc_period(semi_major_axis: f64, grav_param: f64) -> f64 {
    let period: f64 = 2.0 * PI * (semi_major_axis.powi(3)/grav_param).sqrt();
    return period
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
pub fn calc_ecc_vec(
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
/// Inputs
/// ------
/// semi_major_axis: `f64`
///     Semi major axis of orbit
/// 
/// grav_param: `f64`
///     Gravitational parameter of central body
/// 
/// Outputs
/// -------
/// mean_motion: `f64`
///     Mean motion of orbit
pub fn calc_mean_motion(semi_major_axis: f64, grav_param: f64) -> f64 {
    let mean_motion: f64 = 1.0 / (2.0 * PI * (semi_major_axis.powi(3)/grav_param).sqrt());
    mean_motion
}

/// Calculate semi major axis of orbit
/// 
/// Inputs
/// ------
/// grav_param: `f64`
///     Gravitational parameter of central body
/// 
/// mean_motion: `f64`
///     Mean motion of orbit
/// 
/// Outputs
/// -------
/// semi_major_axis: `f64`
///     Semi major axis of orbit
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
/// delta_v_1, delta_v_2 : `f64`
///     Total delta v for maneuver
pub fn calc_hohmann_transfer(
    radius_1: f64,
    radius_2: f64,
    vel_0: f64
) -> (f64, f64) {
    let delta_v_1: f64 = vel_0 * (((2.0 * radius_2)/(radius_1 + radius_2)).sqrt()- 1.0);
    let delta_v_2_num: f64 = (radius_1 / radius_2).sqrt() + (2.0 * radius_1);
    let delta_v_2_den: f64 = (radius_2 *(1.0 + radius_2 / radius_1)).sqrt();
    let delta_v_2: f64 = vel_0 * delta_v_2_num / delta_v_2_den; 
    return (delta_v_1, delta_v_2);
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
        let deltav: (f64, f64) = calc_hohmann_transfer(radius_1, radius_2, vel_0);
        // TODO-TD: fix
        assert_eq!(deltav, (0.23532631438317964, 930.5729285869206));

    }

    #[test]
    fn test_hill_sphere(){
        let earth_mass: f64 = cst::EARTH::MASS;
        let sun_mass: f64 = cst::SUN::MASS;
        let earth_orbit_semi_major: f64 = cst::EARTH::SEMI_MAJOR;
        let earth_orbit_ecc: f64 = cst::EARTH::ORBIT_ECC;

        let sphere_rad: f64 = calc_hill_sphere(
            sun_mass, 
            earth_mass,
            earth_orbit_semi_major,
            earth_orbit_ecc
        );

        assert_eq!(sphere_rad, 7069282679.67562)
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
    
    #[test]
    fn test_calc_lagrange_points()
    {
        let mass_1: f64 = 1.;
        let mass_2: f64 = 0.01;
        let l_points: Vec<Vector3<f64>> = calc_lagrange_points(mass_1, mass_2);

        assert_eq!(l_points[0], Vector3::new(0.8506198417814278, 0., 0.));
        assert_eq!(l_points[1], Vector3::new(1.149380158218572, 0., 0.));
        assert_eq!(l_points[2], Vector3::new(-0.005833333333333334, 0., 0.));
        assert_eq!(l_points[3], Vector3::new(-0.4900990099009901, -0.8660254037844386, 0.));
        assert_eq!(l_points[4], Vector3::new(-0.4900990099009901, 0.8660254037844386, 0.));
    }

    #[test]
    fn check_solar_sytem_params()
    {
        let grav_param: f64 = cst::SUN::GRAV_PARAM;
        let periods: [f64; 8] = [
            0.241,      // Mercury
            0.615,      // Venus
            1.,         // Earth
            1.881,      // Mars
            11.86,      // Jupiter
            29.46,      // Saturn
            84.01,      // Uranus
            164.79];    // Neptune
        
        let semi_majors: [f64; 8] = periods.map(|period| {
            calc_semi_major_axis(grav_param, 1.0 / (period * 365.25 * 86400.))
            }
        );
        assert_eq!(semi_majors, [
            1.9726853165415475e11,  // Mercury Orbit Semi Major Axis
            3.683804104432092e11,   // Venus
            5.09385349788403e11,    // Earth
            7.761976328545394e11,   // Mars
            2.6491277965803633e12,  // Jupiter
            4.858866259527982e12,   // Saturn
            9.770882110650564e12,   // Uranus
            1.5310887052726793e13   // Neptune
        ]);
        
    }
}
