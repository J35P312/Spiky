import sys
import pysam
from collections import defaultdict


# CIGAR operations consuming reference and contributing to coverage
REF_CONSUMING = {0, 7, 8}  # M, =, X


def run(
    bam_path,
    bin_size=1000,
    mapq_threshold=20,
    out=sys.stdout
):
    bam = pysam.AlignmentFile(bam_path, "rb")

    # window_stats[chrom][bin] = [highq_bases, lowq_reads, total_reads]
    window_stats = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        chrom = bam.get_reference_name(read.reference_id)

        is_highq = read.mapping_quality >= mapq_threshold
        is_lowq = not is_highq

        # Track read counts (once per read)
        ref_start = read.reference_start
        ref_end = read.reference_end
        if ref_start is None or ref_end is None:
            continue

        bin_start = ref_start // bin_size
        bin_end = (ref_end - 1) // bin_size

        for b in range(bin_start, bin_end + 1):
            stats = window_stats[chrom][b]
            stats[2] += 1
            if is_lowq:
                stats[1] += 1

        # Only high-MAPQ reads contribute to coverage
        if not is_highq:
            continue

        ref_pos = read.reference_start

        for op, length in read.cigartuples:
            if op in REF_CONSUMING:
                block_start = ref_pos
                block_end = ref_pos + length

                b0 = block_start // bin_size
                b1 = (block_end - 1) // bin_size

                for b in range(b0, b1 + 1):
                    bin_start_bp = b * bin_size
                    bin_end_bp = bin_start_bp + bin_size

                    overlap = max(
                        0,
                        min(block_end, bin_end_bp) -
                        max(block_start, bin_start_bp)
                    )

                    if overlap > 0:
                        window_stats[chrom][b][0] += overlap

                ref_pos += length

            elif op in (2, 3):  # D or N: consume reference, no coverage
                ref_pos += length

            else:
                # I, S, H, P: do not consume reference
                continue

    # Compute global mean coverage (high-MAPQ only)
    total_mean_cov = 0.0
    total_bins = 0

    for chrom in window_stats:
        for b in window_stats[chrom]:
            highq_bases = window_stats[chrom][b][0]
            total_mean_cov += highq_bases / bin_size
            total_bins += 1

    global_mean_cov = total_mean_cov / total_bins if total_bins else 0.0

    # Output BED
    for chrom in sorted(window_stats):
        for b in sorted(window_stats[chrom]):
            start = b * bin_size
            end = start + bin_size

            highq_bases, lowq_reads, total_reads = window_stats[chrom][b]

            mean_cov = highq_bases / bin_size
            frac_cov = mean_cov / global_mean_cov if global_mean_cov > 0 else 0.0
            lowq_frac = lowq_reads / total_reads if total_reads > 0 else 0.0

            out.write(
                f"{chrom}\t{start}\t{end}\t"
                f"{mean_cov:.6f}\t{frac_cov:.6f}\t{lowq_frac:.6f}\n"
            )

    bam.close()
