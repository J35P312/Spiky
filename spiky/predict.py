import sys
import json
import statistics


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

    slope = float(model["slope"])
    intercept = float(model["intercept"])

    # --------------------------------------------------
    # Load regions
    # --------------------------------------------------
    regions = set()
    with open(regions_bed) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            if chrom not in ("Y", "chrY"):
                continue
            regions.add((chrom, int(start), int(end)))

    if not regions:
        raise ValueError("No chrY regions found in regions BED")

    # --------------------------------------------------
    # Extract chrY frac_cov values
    # --------------------------------------------------
    frac_covs = []

    with open(coverage_bed) as bed:
        for line in bed:
            if not line.strip():
                continue
            chrom, start, end, mean_cov, frac_cov, lowq_frac = line.rstrip().split("\t")

            key = (chrom, int(start), int(end))
            if key in regions:
                frac_covs.append(float(frac_cov))

    if not frac_covs:
        raise ValueError("No matching chrY regions found in coverage BED")

    chrY_median_fraction = statistics.median(frac_covs)

    # --------------------------------------------------
    # Predict FFY
    # --------------------------------------------------
    ffy = slope * chrY_median_fraction + intercept

    # --------------------------------------------------
    # Output
    # --------------------------------------------------
    out.write(
        "bed_file\tchrY_median_fraction\tFFY\n"
    )
    out.write(
        f"{coverage_bed}\t{chrY_median_fraction:.6f}\t{ffy:.6f}\n"
    )

