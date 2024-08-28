use anyhow::{Context, Result};
use d4::{
    task::{Task, TaskPartition},
    D4MatrixReader, D4TrackReader, MultiTrackReader,
};
use log::trace;
use std::{panic, path::Path}; // Make sure to include this import for tracing
pub struct TaskPart {
    parent_region: (String, u32, u32),
    thresholds: (f64, f64),
    count_sum: Vec<(u32, u32)>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CallableRegion {
    count: u32,
    begin: u32,
    end: u32,
}

impl<T: Iterator<Item = i32> + ExactSizeIterator> TaskPartition<T> for TaskPart {
    type ParentType = TaskParent;
    type ResultType = Vec<(u32, u32)>;

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
            count_sum: vec![(0, 0); (part_end - part_begin) as usize],
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
            trace!("Value: {:?}", v);
            self.count_sum[left as usize].1 += v as u32;
            trace!(
                "Testing if {:?} is greater than {} and less than {}",
                v,
                self.thresholds.0,
                self.thresholds.1
            );
            if v as f64 >= self.thresholds.0 && v as f64 <= self.thresholds.1 {
                self.count_sum[left as usize].0 += 1
            }
        }
        true
    }

    fn result(&mut self) -> Self::ResultType {
        trace!("Count sum: {:?}", self.count_sum);
        std::mem::take(&mut self.count_sum)
    }
}

pub struct TaskParent {
    chrom: String,
    begin: u32,
    end: u32,
    thresholds: (f64, f64),
    mean_depth_min: f64,
    depth_proportion: f64,
    num_tracks: usize,
    output_counts: bool,
}

impl<T: Iterator<Item = i32> + ExactSizeIterator> Task<T> for TaskParent {
    type Partition = TaskPart;
    type Output = Vec<CallableRegion>;

    fn region(&self) -> (&str, u32, u32) {
        (self.chrom.as_str(), self.begin, self.end)
    }

    fn combine(&self, parts: &[Vec<(u32, u32)>]) -> Vec<CallableRegion> {
        trace!("Combining results from TaskParts: {:?}", parts);
        let mut res = Vec::new();
        let mut current_region: Option<CallableRegion> = None;

        for (index, (count, sum)) in parts.iter().flat_map(|v| v.iter()).enumerate() {
            trace!("count: {}", count);
            let mean = *sum as f64 / self.num_tracks as f64; // mean calculated based on number of individuals total
            let prop = *count as f64 / self.num_tracks as f64; // prop calculated based on number of individuals total
            if mean >= self.mean_depth_min && prop >= self.depth_proportion {
                if let Some(ref mut region) = current_region {
                    // check if current_region is not none
                    // we have a current region, lets see if we can extend it
                    if region.end == (self.begin as usize + index) as u32
                        && (!self.output_counts || region.count == *count)
                    {
                        region.end += 1;
                    } else {
                        // cant extend
                        res.push(region.clone());
                        *region = CallableRegion {
                            count: *count,
                            begin: (self.begin as usize + index) as u32,
                            end: (self.begin as usize + index + 1) as u32,
                        };
                    }
                } else {
                    // no current region lets make one
                    current_region = Some(CallableRegion {
                        count: *count,
                        begin: (self.begin as usize + index) as u32,
                        end: (self.begin as usize + index + 1) as u32,
                    });
                }
            } else if let Some(region) = current_region.take() {
                res.push(region);
            }
        }

        // Push the last region if it exists
        if let Some(region) = current_region {
            res.push(region);
        }

        trace!("Combined result: {:?}", res);
        res
    }
}

pub fn print_as_bed(chrom: &str, begin: u32, output_counts: bool, regions: Vec<CallableRegion>) {
    for region in regions.into_iter() {
        if !output_counts {
            println!(
                "{}\t{}\t{}",
                chrom,
                region.begin + begin,
                region.end + begin,
                
            );
        } else {
            println!(
                "{}\t{}\t{}\t{}",
                chrom,
                region.begin + begin,
                region.end + begin,
                region.count
            );
        }
        
    }
}

pub fn run_d4_tasks(
    d4_file_path: &str,
    thresholds: (f64, f64),
    mean_depth_min: f64,
    depth_proportion: f64,
    output_counts: bool,
) -> Result<()> {
    trace!("Running D4 tasks on file: {}", d4_file_path);

    let tracks: Vec<D4TrackReader> = D4TrackReader::open_tracks(Path::new(d4_file_path), |_| true)
        .context("Failed to open D4 file")?;
    trace!("Successfully opened D4 file with {} tracks", tracks.len());
    let num_tracks = tracks.len();
    let mut reader = D4MatrixReader::new(tracks).unwrap();
    let chrom_regions = reader.chrom_regions();

    trace!("Chromosome regions found: {:?}", chrom_regions);

    let mut tasks = vec![];

    for (chr, begin, end) in chrom_regions.iter() {
        tasks.push(TaskParent {
            chrom: chr.to_string(),
            begin: *begin,
            end: *end,
            thresholds,
            mean_depth_min,
            depth_proportion,
            num_tracks: num_tracks,
            output_counts,
        });
    }

    trace!("Running tasks on D4 reader");
    let result = reader.run_tasks(tasks).unwrap();
    trace!("Tasks completed successfully");

    for r in result.into_iter() {
        print_as_bed(r.chrom, r.begin, output_counts, r.output.to_vec())
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    

    use super::*;

    type RowType = std::vec::IntoIter<i32>;
    fn create_task_parent(output_counts: bool) -> TaskParent {
        TaskParent {
            chrom: String::from("chr1"),
            begin: 0,
            end: 10,
            num_tracks: 10,
            mean_depth_min: 10.0,
            depth_proportion: 0.5,
            thresholds: (0.0, 100.0), // this doesnt matter for combine tests
            output_counts: output_counts,
        }
    }

    #[test]
    fn combine_output_counts_false_adjacent_regions_different_counts() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![vec![(5, 100), (5, 100), (10, 100), (10, 100)]];

        let expected = vec![
            CallableRegion {
                count: 5,
                begin: 0,
                end: 4,
            },
        ];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_output_counts_true_adjacent_regions_different_counts() {
        let output_counts = true;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![vec![(5, 100), (5, 100), (10, 100), (10, 100)]];

        let expected = vec![
            CallableRegion {
                count: 5,
                begin: 0,
                end: 2,
            },
            CallableRegion {
                count: 10,
                begin: 2,
                end: 4,
            },
        ];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }


    #[test]
    fn combine_good_mean_bad_prop() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![vec![(2, 100), (2, 100), (10, 100), (10, 100)]];

        let expected = vec![CallableRegion {
            count: 10,
            begin: 2,
            end: 4,
        }];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_single_continuous_region() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![
            vec![(5, 50), (5, 50), (5, 50)], // mean depth is 5 not callable
        ];

        let expected = vec![];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_non_callable_regions() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![
            vec![(1, 100), (1, 100), (5, 50), (5, 50)], // 0-2 okay mean but bad proportion. 2-4 okay proportion but bad mean
        ];

        let expected = vec![]; // No region should be returned

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_multiple_separated_regions() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![
            vec![(10, 100), (10, 100)], // callable
            vec![(5, 50), (5, 50)],     // not callable
            vec![(0, 0), (0, 0)],       // not callable
            vec![(10, 100), (10, 100)], // callable
        ];

        let expected = vec![
            CallableRegion {
                count: 10,
                begin: 0,
                end: 2,
            },
            CallableRegion {
                count: 10,
                begin: 6,
                end: 8,
            },
        ];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }
}
