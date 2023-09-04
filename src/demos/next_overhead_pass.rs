/*
Calculate next overhead pass windows
*/

use clap::Parser;
use nalgebra::Vector3;
use chrono::{DateTime, Utc};

use donnager::donnager::{
    gravity::kepler as kepler,
    constants as cst, 
    spacetime as xyzt
};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// tle string
    #[arg(short, long)]
    tle_string: String,

    /// Ground station lat lon height
    /// [lat, lon, height]
    #[arg(short, long)]
    ground_station: Vector3<f64>,

    /// Start date of window
    /// [yyyy, mm, dd, hh, mm, ss]
    #[arg(short, long)]
    start_date: Vector6<u32>,

    /// stop date of window
    /// [yyyy, mm, dd, hh, mm, ss]
    #[arg(short, long)]
    stop_date: Vector6<u32>,

}

fn main() {
    
    let args = Args::parse();

    // Orbiting Body
    let tle_str: &str;
    if args.tle_string == None {

        tle_str = "ISS
        1 25544U 98067A   23060.72453421  .00027779  00000-0  50068-3 0  9993
        2 25544  51.6419 140.8390 0005926  43.3718 100.9839 15.49547192385127";
    } else {
        tle_str = args.tle_string
    }
    
    let orbit: kepler::Orbit = kepler::Orbit::from_tle(
        tle_str.to_string()
    );

    let frame: xyzt::ReferenceFrames = xyzt::ReferenceFrames::Cartesian;
    
    // Start time
    let start_date = args.start_date;
    let start_year: i32 = start_date.0 as i32;
    let start_month: u32 = start_date.1;
    let start_day: u32 = start_date.2;
    let start_hour: u32 = start_date.3;
    let start_min: u32 = start_date.4;
    let start_sec: u32 = start_date.5;
    let start_date_time: DateTime<Utc> = xyzt::ymd_hms_to_datetime(
        start_year,
        start_month,
        start_day,
        start_hour,
        start_min,
        start_sec
    );

    // Stop time
    let stop_date = args.stop_date;
    let stop_year: i32 = stop_date.0 as i32;
    let stop_month: u32 = stop_date.1;
    let stop_day: u32 = stop_date.2;
    let stop_hour: u32 = stop_date.3;
    let stop_min: u32 = stop_date.4;
    let stop_sec: u32 = stop_date.5;
    let stop_date_time: DateTime<Utc> = xyzt::ymd_hms_to_datetime(
        stop_year,
        stop_month,
        stop_day,
        stop_hour,
        stop_min,
        stop_sec
    );
    
    // Ground station
    let ground_station: Vector3<f64> = args.ground_station;

    // Calculate ENU

    // TODO-TD: Parallelize
    for date_time in start_date_time..stop_date_time{

        let time_since_epoch: f64 = date_time - orbit.epoch;
        let p_ecef: Vector3<f64> = orbit.calc_motion(
            time_since_epoch, 
            frame, 
            0
        );

        // Calculate angles in ENU


        // Return true each time above horizon
    }


}