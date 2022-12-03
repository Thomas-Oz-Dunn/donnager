mod constants;

fn main() {
    let mass_0: f64 = 5.098e35; // kg
    let mass_1: f64 = 2.4; // kg

    let radius_0: f64 = 344444.0; //m
    let radius_f: f64 = 34444400.0; //m

    let engine_isp: f64 = 543.123; // s

    let delta_v = delta_v(mass_0, radius_0, radius_f);
    let g_0: f64 = constants::GRAV_CONST * mass_0 / radius_0.powf(2.0);
    let mass_fuel = mass_fuel(delta_v, mass_1, engine_isp, g_0);    
    println!("{} kg of fuel", mass_fuel)
}

fn orbital_velocity(
    mass_0: f64,
    radius: f64
) -> f64{
    // KE = PE
    // 1/2 m_1 * v^2 = m_1 * g * r
    // g = F_g / m_1
    // F_g = G * m_1 * m_0 / r^2
    // g = G * m_0 / r^2
    // v = sqrt(2 * g * r)
    // v - sqrt(2 * G * m_0 / r)
    let energy = 2.0 * GRAV_CONST * mass_0 / radius;
    let velocity = energy.sqrt();
    velocity
}

// 1 dimensional, what about hohman transfer and similar?
fn delta_v(
    mass_0: f64,
    radius_0: f64,
    radius_f: f64
) -> f64{
    /*
    mass_0 : float64
        Mass of central body

    radius_0 : float64
        Initial radius
    
    radius_f : float64
        Final velocity
    */
    
    let v_0: f64 = orbital_velocity(mass_0, radius_0);
    let v_f: f64 = orbital_velocity(mass_0, radius_f);
    let delta_v: f64 = v_f - v_0;
    delta_v
}

fn mass_fuel(
    delta_v: f64,
    mass_1: f64,
    engine_isp: f64,
    g_0: f64
) -> f64{
    let power = delta_v / (engine_isp * g_0);
    let mass_f = mass_1 * power.exp();
    mass_f
}