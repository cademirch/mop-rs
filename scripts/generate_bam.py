import pysam
import random
import argparse
import numpy as np


def generate_bam(
    output_bam,
    num_chromosomes,
    chromosome_length,
    max_read_depth,
    fixed_depth=None,
    map_quality=60,
    base_quality=30,
):
    # Cap the read length to 150bp
    read_length = min(chromosome_length, 150)

    # Create a header for the BAM file
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"LN": chromosome_length, "SN": f"chr{i+1}"} for i in range(num_chromosomes)
        ],
    }

    # Open a BAM file for writing
    with pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
        for chrom in range(num_chromosomes):
            chrom_name = f"chr{chrom+1}"

            for pos in range(0, chromosome_length, read_length):
                # Create depth for each chunk (5 bases)
                chunk_depths = []
                for chunk_start in range(0, read_length, 5):
                    # Use fixed depth if specified, otherwise randomize depth for each chunk
                    current_depth = (
                        fixed_depth
                        if fixed_depth is not None
                        else random.randint(0, max_read_depth)
                    )
                    chunk_depths.append(current_depth)

                for chunk_start, current_depth in zip(
                    range(0, read_length, 5), chunk_depths
                ):
                    for _ in range(current_depth):
                        a = pysam.AlignedSegment()
                        a.query_name = f"read_{chrom}_{pos}_{chunk_start}_{random.randint(1, 1000)}"
                        a.query_sequence = "A" * read_length
                        a.flag = 0
                        a.reference_id = chrom
                        a.reference_start = pos
                        a.mapping_quality = map_quality
                        a.cigar = [(0, read_length)]  # Match of read_length bases
                        a.next_reference_id = -1
                        a.next_reference_start = -1
                        a.template_length = read_length
                        a.query_qualities = pysam.qualitystring_to_array(
                            chr(base_quality + 33) * read_length
                        )
                        a.tags = (("NM", 1), ("AS", read_length), ("XS", 0))

                        out_bam.write(a)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate BAM files with controlled depth and quality for testing."
    )

    parser.add_argument("output_bam", type=str, help="Output BAM file path")
    parser.add_argument(
        "-n",
        "--num_chromosomes",
        type=int,
        default=2,
        help="Number of chromosomes to simulate (default: 2)",
    )
    parser.add_argument(
        "-l",
        "--chromosome_length",
        type=int,
        default=1000,
        help="Length of each chromosome in base pairs (default: 1,000)",
    )
    parser.add_argument(
        "-d",
        "--max_read_depth",
        type=int,
        default=30,
        help="Maximum read depth per position (default: 30) when random depth is used",
    )
    parser.add_argument(
        "-D",
        "--fixed_depth",
        type=int,
        help="Fixed read depth per position. If provided, overrides random depth.",
    )
    parser.add_argument(
        "-Q",
        "--map_quality",
        type=int,
        default=60,
        help="Mapping quality of reads (default: 60)",
    )
    parser.add_argument(
        "-q",
        "--base_quality",
        type=int,
        default=30,
        help="Base quality of reads (default: 30)",
    )

    args = parser.parse_args()

    if args.chromosome_length < 10:
        raise ValueError("Chromosome length must be >= 10.")

    generate_bam(
        args.output_bam,
        args.num_chromosomes,
        args.chromosome_length,
        args.max_read_depth,
        args.fixed_depth,
        args.map_quality,
        args.base_quality,
    )
