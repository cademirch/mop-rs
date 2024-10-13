mod d4_utils;
use anyhow::Result;
use clap::{ArgGroup, Parser};
use core::f64;
use d4_utils::run_d4_tasks;
use rayon::ThreadPoolBuilder;

#[derive(Eq, PartialEq, Hash)]
pub struct Region {
    chrom: String,
    begin: u32,
    end: u32,
}

/// A simple tool to process D4 files
#[derive(Parser, Debug)]
#[command(author, version, about="Calculates callable sites from depth statistics.", long_about = None)]
#[command(group(
    ArgGroup::new("input")
        .required(true)
        .args(&["d4_file"]),
))]
struct Args {
    /// The path to the D4 file
    #[arg(long = "d4")]
    d4_file: String,

    /// Minimum depth to consider site callable per individual
    #[arg(short = 'm', long = "min-depth", default_value_t = 0.0)]
    min_depth: f64,

    /// Maximum depth to consider site callable per individual
    #[arg(short = 'M', long = "max-depth", default_value_t = f64::INFINITY)]
    max_depth: f64,

    /// Proporition of samples passing thresholds at site to consider callable
    #[arg(short = 'd', long = "depth-proportion", default_value_t = 1.0)]
    depth_proportion: f64,

    /// Minimum mean depth across all samples at site to consider callable
    #[arg(short = 'u', long = "min-mean-depth", default_value_t = 0.0)]
    mean_depth_min: f64,

    /// Maximum mean depth across all samples at site to consider callable
    #[arg(short = 'U', long = "max-mean-depth", default_value_t = f64::INFINITY)]
    mean_depth_max: f64,
    /// Output number of individuals callable at site. EXPERIMENTAL v0.1.0
    #[arg(short = 'c', long = "output-counts", default_value_t = false)]
    output_counts: bool,

    /// Number of threads to use
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    threads: usize,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    let thresholds = (args.min_depth, args.max_depth);
    let mean_thresholds = (args.mean_depth_min, args.mean_depth_max);
    // Run the D4 tasks with the provided file path and thresholds
    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    run_d4_tasks(
        &args.d4_file,
        thresholds,
        mean_thresholds,
        args.depth_proportion,
        args.output_counts,
    )?;

    Ok(())
}
