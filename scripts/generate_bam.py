import pysam
import random
import argparse


def generate_bam(output_file, n=1, l=10, min_depth=10, max_depth=50):
    assert (
        l >= 10 and l % 10 == 0
    ), "Chromosome length must be at least 10 and divisible by 10."

    chunk_size = l // 5  # Each chunk will be 1/5 of the chromosome length
    reads = []

    # Generate reads for each chromosome
    for chrom_idx in range(n):
        chrom_name = f"chr{chrom_idx + 1}"
        chunks = list(range(0, l, chunk_size))
        random.shuffle(
            chunks
        )  # Randomly shuffle chunks to determine which ones will have depth

        for i, start in enumerate(chunks):
            if i < 3:  # Select 3 out of 5 chunks to have reads
                target_depth = random.randint(min_depth, max_depth)

                for _ in range(target_depth):
                    # Create a synthetic read that spans the entire chunk
                    a = pysam.AlignedSegment()
                    a.query_name = (
                        f"read_{chrom_name}_{start}_{random.randint(1, 1000000)}"
                    )
                    a.query_sequence = (
                        "A" * chunk_size
                    )  # Synthetic sequence matching chunk length
                    a.flag = 0
                    a.reference_id = chrom_idx
                    a.reference_start = start
                    a.mapping_quality = 60
                    a.cigar = [(0, chunk_size)]  # CIGAR: {chunk_size}M
                    a.next_reference_id = -1
                    a.next_reference_start = -1
                    a.template_length = 0
                    a.query_qualities = pysam.qualitystring_to_array("I" * chunk_size)

                    # Add the read to the list
                    reads.append(a)

    # Sort reads by reference_id and reference_start
    reads.sort(key=lambda x: (x.reference_id, x.reference_start))

    # Open a new BAM file for writing
    with pysam.AlignmentFile(
        output_file,
        "wb",
        header={
            "HD": {"VN": "1.0"},
            "SQ": [{"LN": l, "SN": f"chr{i+1}"} for i in range(n)],
        },
    ) as bam_file:
        # Write the sorted reads to the BAM file
        for read in reads:
            bam_file.write(read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a BAM file with specified parameters."
    )
    parser.add_argument("output_file", help="The output BAM file name.")
    parser.add_argument(
        "-n", type=int, default=1, help="Number of chromosomes (default: 1)."
    )
    parser.add_argument(
        "-l", type=int, default=10, help="Length of each chromosome (default: 10)."
    )
    parser.add_argument(
        "--min_depth",
        type=int,
        default=10,
        help="Minimum depth of coverage (default: 10).",
    )
    parser.add_argument(
        "--max_depth",
        type=int,
        default=50,
        help="Maximum depth of coverage (default: 50).",
    )

    args = parser.parse_args()

    generate_bam(
        args.output_file,
        n=args.n,
        l=args.l,
        min_depth=args.min_depth,
        max_depth=args.max_depth,
    )
