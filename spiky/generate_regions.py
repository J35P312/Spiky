import sys
import statistics
from collections import defaultdict
import numpy as np


def run(
    sample_list_path,
    female_median_threshold=0.01,
    out_bed=sys.stdout,
    out_csv_path=None
):
    """
    Generate robust chrY regions for downstream analysis.
    """

    if out_csv_path is None:
        raise ValueError("CSV output path must be provided")

    # Load sample list
    samples = []
    with open(sample_list_path) as f:
        for line in f:
            if not line.strip():
                continue
            bed_path, sex = line.rstrip().split(",")
            samples.append((bed_path, sex))

    # Data structures
    # bins[(chrom, start, end)][bed_path] = frac_cov
    bins = defaultdict(dict)

    # Read coverage BEDs
    for bed_path, sex in samples:
        with open(bed_path) as bed:
            for line in bed:
                if not line.strip():
                    continue

                chrom, start, end, mean_cov, frac_cov, lowq_frac = line.rstrip().split("\t")

                if chrom not in ("Y", "chrY"):
                    continue

                key = (chrom, int(start), int(end))
                bins[key][bed_path] = float(frac_cov)

    # Prepare CSV output
    csv_out = open(out_csv_path, "w")
    csv_out.write("position,bed_file,sex,status\n")

    # Process bins
    for (chrom, start, end), values in sorted(bins.items()):
        male_vals = []
        female_vals = []

        for bed_path, sex in samples:
            if bed_path not in values:
                continue
            frac_cov = values[bed_path]
            if sex == "male":
                male_vals.append(frac_cov)
            else:
                female_vals.append(frac_cov)

        # Skip bins without both groups
        if not male_vals or not female_vals:
            continue

        male_vals_np = np.array(male_vals)
        female_vals_np = np.array(female_vals)

        male_p5 = np.percentile(male_vals_np, 5)
        female_p95 = np.percentile(female_vals_np, 95)

        male_median = statistics.median(male_vals)
        female_median = statistics.median(female_vals)

        status = (
            "ok"
            if male_p5 > female_p95 and female_median < female_median_threshold
            else "fail"
        )

        # Write BED (one line per bin)
        if status == "ok":
            out_bed.write(
                f"{chrom}\t{start}\t{end}\t"
                f"{male_median:.6f}\t{female_median:.6f}\t"
                f"{male_p5:.6f}\t{female_p95:.6f}\n"
            )
      
        # Write CSV (one line per sample per bin)
        for bed_path, sex in samples:
            csv_out.write(
                f"{start},{bed_path},{sex},{status}\n"
            )

    csv_out.close()
