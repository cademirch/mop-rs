mod d4_utils;

use clap::Parser;
use d4_utils::run_d4_tasks;
use anyhow::Result;

/// A simple tool to process D4 files
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// The path to the D4 file
    #[arg(short, long)]
    d4_file: String,

    // /// Threshold values
    // #[arg(short, long, default_value = (0.5, 1.0))]
    // thresholds: (f64, f64),
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    let thresholds = (1.0, 10.0);
    let d = 1.0;
    let y = 0.0;
    // Run the D4 tasks with the provided file path and thresholds
    run_d4_tasks(&args.d4_file, thresholds, d ,y )?;

    println!("Processing complete!");

    Ok(())
}
