import random
import pysam
import os
import numpy as np


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
    output_dir="bam_simulation",
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
    for bam_index in range(num_bams):
        bam_file_path = os.path.join(output_dir, f"simulated_{bam_index + 1}.bam")

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
    return truth_depth_sum, truth_depth_count, intervals


def make_truth_file(
    truth_sum,
    truth_coumt,
    num_chromosomes,
    min_mean_depth,
    proportion,
    num_bams,
    intervals,
):
    for chr_num in range(num_chromosomes):
        chromosome = f"chr{chr_num + 1}"
        sums = truth_sum[chromosome]
        counts = truth_count[chromosome]

        for i, (sum, count) in enumerate(zip(sums, counts)):
            mean = sum / num_bams
            prop = count / num_bams
            start, stop = intervals[i]
            if mean >= min_mean_depth and prop >= proportion:
                print(f"{chromosome}\t{start}\t{stop}")


if __name__ == "__main__":
    num_bams = 3
    num_chromosomes = 2
    chr_length = 1000
    min_depth = 10
    max_depth = 30
    prop_outside_min_max = 0.3
    min_mean_depth = 1
    proportion = 1

    truth_sum, truth_count, intervals = generate_bam_files(
        num_bams=num_bams,
        num_chromosomes=num_chromosomes,
        chr_length=chr_length,
        min_depth=min_depth,
        max_depth=max_depth,
        interval_mean_size=100,
        interval_size_std=20,
    )

    make_truth_file(
        truth_sum,
        truth_count,
        num_chromosomes,
        min_mean_depth,
        proportion,
        num_bams,
        intervals,
    )

#     print("BAM files and corresponding merged BED files have been generated.")
