use anyhow::{Context, Result};
use d4::{
    task::{Mean, Task, TaskOutputVec, TaskPartition},
    D4MatrixReader, D4TrackReader, MultiTrackReader,
};
use log::{error, trace};
use std::{panic, path::Path}; // Make sure to include this import for tracing

pub struct TaskPart {
    parent_region: (String, u32, u32),
    part_begin: u32,
    part_end: u32,
    thresholds: (f64, f64),
    count_sum: Vec<(u32,u32)>,
    mean_depth_min: f64,
    depth_proportion: f64
    
}

impl<T: Iterator<Item = i32> + ExactSizeIterator> TaskPartition<T> for TaskPart {
    type ParentType = TaskParent;
    type ResultType = Vec<(u32,u32)>;

    fn new(part_begin: u32, part_end: u32, parent: &Self::ParentType) -> Self {
        trace!(
            "Creating new TaskPart for region: {}-{}, in parent: {}:{}-{}",
            part_begin,
            part_end,
            parent.chrom,
            parent.begin,
            parent.end
        );
        Self {
            parent_region: (parent.chrom.clone(), parent.begin, parent.end),
            thresholds: parent.thresholds,
            part_begin: part_begin,
            part_end: part_end,
            count_sum: vec![(0,0); (part_end - part_begin) as usize],
            mean_depth_min: parent.mean_depth_min,
            depth_proportion: parent.depth_proportion
            
        }
    }

    fn feed_range(&mut self, left: u32, right: u32, value: &mut T) -> bool {
        trace!(
            "TaskPart(parent: {:?}) Feeding range: {}-{}",
            self.parent_region,
            left,
            right
        );

        if right - left > 1 {
            panic!("Invalid range: {}-{} in task partition in parent: {:?}. The difference between right and left must be 1.", left, right, self.parent_region);
        }

        for v in value.into_iter() {
            self.count_sum[left as usize].1 += v as u32;

            if self.mean_depth_min > 0 as f64 {
                self.count_sum[left as usize].0 += 1
            }
            else if v as f64 >= self.thresholds.0 && v as f64 <= self.thresholds.1 {
                self.count_sum[left as usize].0 += 1
            }
        }
        true
    }

    fn result(&mut self) -> Self::ResultType {
        std::mem::take(& mut self.count_sum)
    }
}

pub struct TaskParent {
    chrom: String,
    begin: u32,
    end: u32,
    thresholds: (f64, f64),
    mean_depth_min: f64,
    depth_proportion: f64,
    num_tracks: usize
    
}


impl<T: Iterator<Item = i32> + ExactSizeIterator> Task<T> for TaskParent {
    type Partition = TaskPart;
    type Output = Vec<(bool, u32)>;

    fn region(&self) -> (&str, u32, u32) {
        (self.chrom.as_str(), self.begin, self.end)
    }

    fn combine(&self, parts: &[Vec<(u32, u32)>]) -> Self::Output {
        trace!("Combining results from TaskParts: {:?}", parts);
        let mut res = vec![(false, 0); self.end as usize - self.begin as usize];
        // mean_depth_min and depth_proporiton need to be mutually exclusive. if neither, then output counts within thresholds.
        for (index, (count, sum)) in parts.iter().flat_map(|v| v.iter()).enumerate() {
            if self.mean_depth_min > 0.0 {
                let mean = *sum as f64 / *count as f64;
                if mean > self.mean_depth_min {
                    res[index].0 = true;
                    trace!("Mean depth: {}", mean)
                }
            } else if self.depth_proportion >= 0.0 {
                let prop = *count as f64 / self.num_tracks as f64;
                if prop >= self.depth_proportion {
                    res[index].0 = true
                }

            } else {
                res[index].1 = *count
            }
        }
       

        trace!("Combined result: {:?}", res);
        res
    }
}

pub fn run_d4_tasks(d4_file_path: &str, thresholds: (f64, f64), mean_depth_min: f64, depth_proportion: f64) -> Result<TaskOutputVec<Vec<(bool, u32)>>> {
    trace!("Running D4 tasks on file: {}", d4_file_path);

    let mut tracks: Vec<D4TrackReader> =
        D4TrackReader::open_tracks(Path::new(d4_file_path), |_| true)
            .context("Failed to open D4 file")?;
    trace!("Successfully opened D4 file with {} tracks", tracks.len());
    let num_tracks = tracks.len();
    let mut reader = D4MatrixReader::new(tracks).unwrap();
    let chrom_regions = reader.chrom_regions();
    trace!("Chromosome regions found: {:?}", chrom_regions);

    let mut tasks = vec![];

    for (chr, begin, end) in chrom_regions.iter() {
        tasks.push(TaskParent {
            chrom: chr.to_string(), // Assuming `chr` is a reference type like `&str` or `String`
            begin: *begin,
            end: *end,
            thresholds,
            mean_depth_min,
            depth_proportion,
            num_tracks: num_tracks
        });
    }

    trace!("Running tasks on D4 reader");
    let result = reader.run_tasks(tasks).unwrap();
    trace!("Tasks completed successfully");

    Ok(result)
}