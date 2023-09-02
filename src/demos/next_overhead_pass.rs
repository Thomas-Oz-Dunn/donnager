/*
Calculate next overhead pass windows
*/

use clap::Parser;

use nalgebra::Vector3;

use donnager::donnager::{
    gravity::kepler as kepler,
    constants as cst, 
    spacetime as xyzt};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// tle string
    #[arg(short, long)]
    tle_string: String,

    /// Ground station lat lon height
    #[arg(short, long)]
    ground_station: Vector3<f64>,

    /// start date of window
    #[arg(short, long)]
    start_date: Vector3<u32>,

    /// stop date of window
    #[arg(short, long)]
    stop_date: Vector3<u32>,


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

    let frame = xyzt::ReferenceFrames::Cartesian;
    
    // Start time

    // Stop time

    
    let p_ecef: Vector3<f64> = orbit.calc_motion(time_since_epoch, frame, 0);
    
    // Ground station
    let ground_station: Vector3<f64> = args.ground_station;

    


}