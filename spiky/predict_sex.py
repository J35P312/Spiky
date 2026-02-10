import sys
import pysam


def run(
    coverage_bed_path,
    reference_fasta,
    y_frac_threshold=0.02,
    lowq_threshold=0.1,
    max_n_frac=0.2,
    out=sys.stdout
):
    """
    Predict sex based on fractional coverage on chromosome Y.
    """

    fasta = pysam.FastaFile(reference_fasta)

    frac_cov_values = []

    with open(coverage_bed_path) as bed:
        for line in bed:
            if not line.strip():
                continue

            chrom, start, end, mean_cov, frac_cov, lowq_frac = line.rstrip().split("\t")

            if chrom not in ("Y", "chrY"):
                continue

            lowq_frac = float(lowq_frac)
            if lowq_frac >= lowq_threshold:
                continue

            start = int(start)
            end = int(end)
            frac_cov = float(frac_cov)

            # Fetch reference sequence and compute N fraction
            seq = fasta.fetch(chrom, start, end).upper()
            if not seq:
                continue

            n_frac = seq.count("N") / len(seq)
            if n_frac >= max_n_frac:
                continue

            frac_cov_values.append(frac_cov)

    fasta.close()

    if frac_cov_values:
        mean_y_frac_cov = sum(frac_cov_values) / len(frac_cov_values)
    else:
        mean_y_frac_cov = 0.0

    prediction = "male" if mean_y_frac_cov >= y_frac_threshold else "female"

    out.write(f"{coverage_bed_path},{prediction}\n")
