use clap::Parser;
use rand::Rng;
use rand_distr::{Distribution, Normal};
use rust_htslib::bam::header::{self, Header, HeaderRecord};
use rust_htslib::bam::record::{Cigar, CigarString, Record};
use rust_htslib::bam::{Format, Writer};
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Number of BAM files to generate
    #[arg(short, long)]
    num_bams: usize,

    /// Number of chromosomes
    #[arg(short, long)]
    num_chrs: usize,

    /// Length of each chromosome
    #[arg(short, long)]
    chr_length: usize,

    /// Minimum depth of coverage
    #[arg(short, long)]
    min_depth: usize,

    /// Maximum depth of coverage
    #[arg(short, long)]
    max_depth: usize,

    /// Minimum mean depth for intervals to be included in the truth file
    #[arg(short, long)]
    min_mean_depth: f64,

    /// Minimum proportion of BAM files covering an interval for it to be included in the truth file
    #[arg(short, long)]
    proportion: f64,

    /// Output directory for BAM files and the truth file
    #[arg(short, long)]
    outdir: String,
}

fn generate_intervals(total_len: usize, mean_size: f32, std_dev: f32) -> Vec<(usize, usize)> {
    let mut intervals = Vec::new();
    let mut current_pos = 0;
    let distr = Normal::new(mean_size, std_dev).unwrap();
    while current_pos < total_len {
        let mut size = distr.sample(&mut rand::thread_rng()) as i32;
        size = max(1, min(size, (total_len - current_pos) as i32));
        let interval = (current_pos, (current_pos + size as usize));
        intervals.push(interval);
        current_pos += size as usize;
    }
    intervals
}

fn write_bam(
    num_bams: usize,
    min_depth: usize,
    max_depth: usize,
    intervals: &Vec<(usize, usize)>,
    num_chrs: usize,
    chr_length: usize,
    outdir: &str,
) -> Result<HashMap<usize, Vec<(usize, usize)>>, Box<dyn std::error::Error>> {
    let mut res = HashMap::new();
    let mut bam_header = Header::new();

    for chr_index in 0..=num_chrs {
        let sum_count = vec![(0, 0); chr_length];
        res.insert(chr_index, sum_count);
        let mut rec = HeaderRecord::new(b"SQ");
        rec.push_tag(b"SN", &format!("sq{}", chr_index).to_string())
            .push_tag(b"LN", &chr_index);
        bam_header.push_record(&rec);
    }

    std::fs::create_dir_all(&outdir)?;
    let bamlist_path = Path::new(&outdir).join("bamlist.txt");
    let mut bamlist = File::create(&bamlist_path)?;

    for bam_index in 0..num_bams {
        let file_path = Path::new(&outdir).join(format!("{}.bam", bam_index));
        writeln!(bamlist, "{}", file_path.display())?;
        let mut writer = Writer::from_path(&file_path, &bam_header, Format::Bam)?;
        let _ = writer.set_threads(8);
        let mut rng = rand::thread_rng();
        for chr_index in 0..num_chrs {
            for (interval_index, (start, stop)) in intervals.iter().enumerate() {
                let good_depth = rng.gen_bool(0.5);
                let mut depth = 0;
                if let Some(sum_count_vec) = res.get_mut(&chr_index) {
                    if let Some(sum_count) = sum_count_vec.get_mut(interval_index) {
                        if good_depth {
                            depth = rng.gen_range(min_depth..=max_depth);
                            sum_count.1 += 1;
                        } else {
                            let lo = rng.gen_bool(0.5);
                            depth = if lo {
                                rng.gen_range(0..=(min_depth - 1))
                            } else {
                                rng.gen_range((max_depth + 1)..(max_depth + 1000))
                            };
                        }
                        sum_count.0 += depth;
                    }
                }
                let cigar = CigarString(vec![Cigar::Match(*stop as u32 - *start as u32)]);
                let sequence = vec![b'A'; stop - start];
                let quality = vec![b'I'; stop - start];

                for i in 0..depth {
                    let mut rec = Record::new();
                    let qname = format!("read_{}_pos_{}", i, start);
                    rec.set(
                        qname.as_bytes(),
                        Some(&cigar),
                        &sequence,
                        &quality,
                    );
                    rec.set_pos(*start as i64);
                    rec.set_tid(0);
                    rec.set_flags(0);
                    writer.write(&rec)?;
                }
            }
        }
    }
    Ok(res)
}

fn make_truth_file(
    truth: &HashMap<usize, Vec<(usize, usize)>>,
    num_chromosomes: usize,
    min_mean_depth: f64,
    proportion: f64,
    num_bams: usize,
    intervals: &[(usize, usize)],
    outdir: &str,
) -> io::Result<()> {
    let truth_path = Path::new(outdir).join("truth.bed");
    let mut file = File::create(truth_path)?;

    for chr_num in 0..num_chromosomes {
        let chromosome = format!("sq{}", chr_num + 1);
        if let Some(data) = truth.get(&chr_num) {
            let mut merged_intervals: Vec<(usize, usize)> = Vec::new();
            let mut current_start: Option<usize> = None;
            let mut current_end: Option<usize> = None;

            for (i, ((sum, count), (start, stop))) in data.iter().zip(intervals.iter()).enumerate()
            {
                let mean = *sum as f64 / num_bams as f64;
                let prop = *count as f64 / num_bams as f64;
                if mean >= min_mean_depth && prop >= proportion {
                    match (current_start, current_end) {
                        (None, None) => {
                            current_start = Some(*start);
                            current_end = Some(*stop);
                        }
                        (Some(_), Some(end)) if *start == end => {
                            current_end = Some(*stop);
                        }
                        (Some(start_val), Some(end_val)) => {
                            merged_intervals.push((start_val, end_val));
                            current_start = Some(*start);
                            current_end = Some(*stop);
                        }
                        _ => {}
                    }
                }
            }

            if let (Some(start), Some(end)) = (current_start, current_end) {
                merged_intervals.push((start, end));
            }

            for (start, stop) in merged_intervals {
                writeln!(file, "{}\t{}\t{}", chromosome, start, stop)?;
            }
        }
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let intervals = generate_intervals(args.chr_length, 100.0, 20.0);
    let truth = write_bam(
        args.num_bams,
        args.min_depth,
        args.max_depth,
        &intervals,
        args.num_chrs,
        args.chr_length,
        &args.outdir,
    )?;
    let _ = make_truth_file(
        &truth,
        args.num_chrs,
        args.min_mean_depth,
        args.proportion,
        args.num_bams,
        &intervals,
        &args.outdir,
    );
    Ok(())
}
