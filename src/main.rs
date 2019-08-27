#![feature(non_ascii_idents)]

use std::time::Instant;

/// the following function enacts the algorithm as specified in the assignment
fn advection(N: usize, dt: f32) -> f32 {
    let now = Instant::now(); // a time instant used for benchmarking

    let tmax = 2.0;
    let nbstep = (tmax / dt) as isize;
    let xmin = 0.0;
    let xmax = 1.0;
    let dx = (xmax - xmin) / 2.0;
    
    let v = 1.0;
    let xc = 0.25;

    let α = v * dt/(2.0 * dx);

    let x: Vec<f32> = (0..=N+2)
                            .map(|i| xmin + (i as f32 - 1.0))
                            .collect();
    
    let u0: Vec<f32> = (0..=N+2)
                    .map(|i| -200.0 * (x[i] - xc).powi(2))
                    .collect();

    let mut u = u0.clone();
    let mut unew = u0.clone();

    for timestep in 1..=nbstep {
        let current_time = timestep as f32 * dt;
        
        for j in 1..=(N + 1) {
            unew[j] = u[j] - α * (u[j+1] - u[j-1] + 0.5 * (u[j+1] - 2.0*u[j] + u[j-1]))
        }
        u = unew.clone();

        u[0] = u[N+1];
        u[N+2] = u[1]
    }

    // measuring and returing elasped time
    now.elapsed().as_secs_f32()
}

/// a main function similar to java's or c's main()
fn main() {
    // runs the calculations 10 tines in a loop and benchmarks for case 1
    let case_1_results = (0..10).map(|_| advection(103, 0.0009))
                                .fold(0.0, |sum, x| sum + x) / 10.0; 

    // runs the calculations 10 tines in a loop and benchmarks for case 2
    let case_2_results = (0..10).map(|_| advection(1003, 0.00009))
                                .fold(0.0, |sum, x| sum + x) / 10.0; 
        
    println!("Average time taken for Case 1: N = 103, dt = 0.0009 over 10 runs is {} seconds", case_1_results);
    println!("Average time taken for Case 2: N = 1003, dt = 0.00009 over 10 runs is {} seconds", case_2_results);
}