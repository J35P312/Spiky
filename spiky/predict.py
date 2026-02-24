import sys
import json
import statistics
import numpy as np

def compute_median_cov_per_sample(coverage_bed_file, ok_bins):
    """
    Compute median coverage for a given chromosome in a sample coverage BED.
    """
    cov_values = []

    with open(coverage_bed_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            chr_ = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            cov = float(parts[4])  # coverage fraction column in coverage BED

            if (chr_, start, end) in ok_bins:
                cov_values.append(cov)

    return np.median(cov_values) if cov_values else 0.0

def load_ok_bins(regions_bed, chromosome):
    """
    Load OK bins from generate_regions BED for a given chromosome.
    """
    ok_bins = set()
    with open(regions_bed) as f:
        for line in f:
            parts = line.strip().split('\t')
            chr_ = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chr_ == chromosome:
                ok_bins.add((chr_, start, end))
    return ok_bins


def run(
    coverage_bed,
    regions_bed,
    model_file,
    out=sys.stdout
):
    """
    Predict fetal fraction from chrY coverage using a trained model.
    """

    # --------------------------------------------------
    # Load model
    # --------------------------------------------------
    with open(model_file) as f:
        model = json.load(f)

    slopeY = model["FFY"]["slope"]
    interceptY = model["FFY"]["intercept"]
    slopeX = model["FFX"]["slope"]
    interceptX = model["FFX"]["intercept"]

    # --------------------------------------------------
    # Load regions
    # --------------------------------------------------
    ok_bins_Y = load_ok_bins(regions_bed, "Y")
    ok_bins_X = load_ok_bins(regions_bed, "X")

    # --------------------------------------------------
    # Extract chrY frac_cov values
    # --------------------------------------------------
    lines = ["Sample\tmedian_chrY\tFFY\tFFX"]

    median_chrY = compute_median_cov_per_sample(coverage_bed, ok_bins_Y)
    median_chrX = compute_median_cov_per_sample(coverage_bed, ok_bins_X)

    FFY_pred = slopeY * median_chrY + interceptY
    FFX_pred = slopeX * median_chrX + interceptX

    lines.append(f"{coverage_bed}\t{median_chrY:.4f}\t{FFY_pred:.4f}\t{FFX_pred:.4f}")

    output_text = "\n".join(lines)
    # --------------------------------------------------
    # Output
    # --------------------------------------------------
    print(output_text)
