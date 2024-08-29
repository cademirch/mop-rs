import random
import pysam
import os
import numpy as np
from snakemake.script import snakemake


def generate_random_intervals(total_length=1000, mean_size=None, std_dev=50):
    if mean_size is None:
        mean_size = total_length / 2

    intervals = []
    current_position = 0

    while current_position < total_length:
        # Generate a random interval size from a normal distribution
        size = int(np.random.normal(mean_size, std_dev))

        # Ensure the size is positive and does not exceed the remaining length
        size = max(1, min(size, total_length - current_position))

        # Create the interval as a tuple (start, end)
        interval = (current_position, current_position + size)

        # Add the interval to the list
        intervals.append(interval)

        # Move the current position to the end of the current interval
        current_position += size

    return intervals


def generate_bam_files(
    num_bams,
    num_chromosomes,
    chr_length,
    min_depth,
    max_depth,
    interval_mean_size,
    interval_size_std,
    output_dir="simulated_bams",
):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    intervals = generate_random_intervals(
        chr_length, interval_mean_size, interval_size_std
    )
    truth_depth_sum = {
        f"chr{chr_num + 1}": [0] * len(intervals) for chr_num in range(num_chromosomes)
    }

    truth_depth_count = {
        f"chr{chr_num + 1}": [0] * len(intervals) for chr_num in range(num_chromosomes)
    }
    with open(os.path.join(output_dir, "bamlist.txt"), "w") as bamlist:
        for bam_index in range(num_bams):
            bam_file_path = os.path.join(output_dir, f"simulated_{bam_index + 1}.bam")
            bamlist.write(bam_file_path + "\n")
            header = {
                "HD": {"VN": "1.0"},
                "SQ": [
                    {"LN": chr_length, "SN": f"chr{chr_num + 1}"}
                    for chr_num in range(num_chromosomes)
                ],
            }
            with pysam.AlignmentFile(bam_file_path, "wb", header=header) as bam_file_obj:
                for chr_num in range(num_chromosomes):
                    chromosome = f"chr{chr_num + 1}"

                    for i, (start, stop) in enumerate(intervals):
                        good_depth = bool(random.getrandbits(1))

                        if good_depth:
                            depth = random.randint(min_depth, max_depth)

                            truth_depth_count[chromosome][i] += 1

                        else:
                            lo = bool(random.getrandbits(1))
                            if lo:
                                depth = random.randint(0, min_depth - 1)
                            else:
                                depth = random.randint(max_depth + 1, max_depth + 1000)
                        truth_depth_sum[chromosome][i] += depth
                        # print(f"{bam_file_path}, {chromosome}, {(start,stop)}, {depth}")

                        for _ in range(depth):

                            a = pysam.AlignedSegment()
                            a.query_name = f"read_{bam_index + 1}_{start + 1}"
                            a.query_sequence = "A" * (stop - start)
                            a.flag = 0
                            a.reference_id = chr_num
                            a.reference_start = start
                            a.mapping_quality = 255
                            a.cigar = [(0, stop - start)]
                            a.next_reference_id = -1
                            a.next_reference_start = -1
                            a.template_length = 0
                            a.query_qualities = pysam.qualitystring_to_array(
                                "I" * (stop - start)
                            )

                            bam_file_obj.write(a)
            pysam.index(bam_file_path)
    return truth_depth_sum, truth_depth_count, intervals


def make_truth_file(
    truth_sum,
    truth_count,
    num_chromosomes,
    min_mean_depth,
    proportion,
    num_bams,
    intervals,
    output_dir="simulated_bams"
):
    with open(os.path.join(output_dir, "truth.bed"), "w") as f:
        for chr_num in range(num_chromosomes):
            chromosome = f"chr{chr_num + 1}"
            sums = truth_sum[chromosome]
            counts = truth_count[chromosome]

            # Initialize merged intervals list
            merged_intervals = []
            
            # Temporary storage for merging intervals
            current_start, current_end = None, None
            
            for i, (sum, count) in enumerate(zip(sums, counts)):
                mean = sum / num_bams
                prop = count / num_bams
                start, stop = intervals[i]

                if mean >= min_mean_depth and prop >= proportion:
                    # Check if we need to merge intervals
                    if current_start is None:
                        # Initialize the first interval
                        current_start, current_end = start, stop
                    elif start == current_end:
                        # Merge the current interval with the previous one
                        current_end = stop
                    else:
                        # Write the previous merged interval and start a new one
                        merged_intervals.append((current_start, current_end))
                        current_start, current_end = start, stop

            # Append the last interval after the loop
            if current_start is not None:
                merged_intervals.append((current_start, current_end))

            # Write the merged intervals to the BED file
            for start, stop in merged_intervals:
                f.write(f"{chromosome}\t{start}\t{stop}\n")


if __name__ == "__main__":
    num_bams = snakemake.params["num_bams"]
    num_chromosomes = snakemake.params["num_chrs"]
    chr_length = snakemake.params["chr_length"]
    min_depth = snakemake.params["min_depth"]
    max_depth = snakemake.params["max_depth"]
    min_mean_depth = snakemake.params["min_mean_depth"]
    proportion = snakemake.params["proportion"]
    outdir = snakemake.params["outdir"]
    truth_sum, truth_count, intervals = generate_bam_files(
        num_bams=num_bams,
        num_chromosomes=num_chromosomes,
        chr_length=chr_length,
        min_depth=min_depth,
        max_depth=max_depth,
        interval_mean_size=100,
        interval_size_std=20,
        output_dir=outdir
    )

    make_truth_file(
        truth_sum,
        truth_count,
        num_chromosomes,
        min_mean_depth,
        proportion,
        num_bams,
        intervals,
        output_dir=outdir
    )

