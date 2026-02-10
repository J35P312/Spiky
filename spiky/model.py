import sys
import json
import statistics
import numpy as np


def run(
    training_csv,
    regions_bed,
    model_out,
    out_csv=sys.stdout
):
    """
    Train linear regression model predicting fetal fraction
    from chrY median fractional coverage.
    """

    # --------------------------------------------------
    # Read training table
    # --------------------------------------------------
    samples = []
    with open(training_csv) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip().split(",")
            if parts[0].lower().endswith(".bed"):
                bed_path = parts[0]
                ff = float(parts[1])
                samples.append((bed_path, ff))

    if not samples:
        raise ValueError("No valid samples found in training CSV")

    # --------------------------------------------------
    # Read regions
    # --------------------------------------------------
    regions = []
    with open(regions_bed) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            if chrom not in ("Y", "chrY"):
                continue
            regions.append((chrom, int(start), int(end)))

    if not regions:
        raise ValueError("No chrY regions found in regions BED")

    # --------------------------------------------------
    # Extract median chrY frac_cov per sample
    # --------------------------------------------------
    x_vals = []  # chrY median frac_cov
    y_vals = []  # fetal fraction

    per_sample_median = {}

    for bed_path, ff in samples:
        frac_covs = []

        with open(bed_path) as bed:
            for line in bed:
                if not line.strip():
                    continue
                chrom, start, end, mean_cov, frac_cov, lowq_frac = line.rstrip().split("\t")

                key = (chrom, int(start), int(end))
                if key in regions:
                    frac_covs.append(float(frac_cov))

        if not frac_covs:
            raise ValueError(f"No matching regions found in {bed_path}")

        median_frac = statistics.median(frac_covs)
        per_sample_median[bed_path] = median_frac

        x_vals.append(median_frac)
        y_vals.append(ff)

    x = np.array(x_vals)
    y = np.array(y_vals)

    # --------------------------------------------------
    # Linear regression: y = slope * x + intercept
    # --------------------------------------------------
    slope, intercept = np.polyfit(x, y, 1)

    # --------------------------------------------------
    # Write model (JSON)
    # --------------------------------------------------
    model = {
        "slope": float(slope),
        "intercept": float(intercept)
    }

    with open(model_out, "w") as f:
        json.dump(model, f, indent=2)

    # --------------------------------------------------
    # Write evaluation CSV
    # --------------------------------------------------
    out_csv.write(
        "bed_file,fetal_fraction,predicted_fetal_fraction,chrY_median_fraction\n"
    )

    for bed_path, ff in samples:
        median_frac = per_sample_median[bed_path]
        predicted = slope * median_frac + intercept

        out_csv.write(
            f"{bed_path},{ff:.6f},{predicted:.6f},{median_frac:.6f}\n"
        )
