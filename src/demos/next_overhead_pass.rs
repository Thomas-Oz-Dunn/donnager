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
    /// Ground station
    #[arg(short, long)]
    ground_station: Vector3<f64>,

}

fn main() {
    
    let args = Args::parse();

    // Orbiting Body
    let tle_str: &str = "ISS
    1 25544U 98067A   23060.72453421  .00027779  00000-0  50068-3 0  9993
    2 25544  51.6419 140.8390 0005926  43.3718 100.9839 15.49547192385127";
    
    let orbit: kepler::Orbit = kepler::Orbit::from_tle(
        tle_str.to_string()
    );

    let frame = xyzt::ReferenceFrames::Cartesian;
    
    // Start time

    // Stop time

    // Ground station
    
    let p_ecef = orbit.calc_motion(time_since_epoch, frame, 0);

    let ground_station: Vector3<f64> = args.ground_station;

    


}