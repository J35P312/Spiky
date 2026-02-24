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
            if chrom not in ("Y", "chrY","X","chrX"):
                continue
            regions.append((chrom, int(start), int(end)))

    if not regions:
        raise ValueError("No chrY regions found in regions BED")

    # --------------------------------------------------
    # Extract median chrY frac_cov per sample
    # --------------------------------------------------
    x_ffy_vals = []  # chrY median frac_cov
    x_ffx_vals = []  # chrY median frac_cov
    y_vals = []  # fetal fraction

    per_sample_median_y = {}
    per_sample_median_x = {}

    for bed_path, ff in samples:
        frac_covs_y = []
        frac_covs_x = []

        with open(bed_path) as bed:
            for line in bed:
                if not line.strip():
                    continue
                chrom, start, end, mean_cov, frac_cov, lowq_frac = line.rstrip().split("\t")
                key = (chrom, int(start), int(end))
                if key in regions:
                    if chrom in ("Y","chrY"):
                        frac_covs_y.append(float(frac_cov))
                    elif chrom in ("X","chrX"):
                        frac_covs_x.append(float(frac_cov))


        median_frac_y = statistics.median(frac_covs_y)
        median_frac_x = statistics.median(frac_covs_x)
   
        per_sample_median_y[bed_path] = median_frac_y
        per_sample_median_x[bed_path] = median_frac_x

        x_ffy_vals.append(median_frac_y)
        x_ffx_vals.append(median_frac_x)
        y_vals.append(ff)

    x_ffy = np.array(x_ffy_vals)
    x_ffx = np.array(x_ffx_vals)
    y = np.array(y_vals)

    # Fit linear regression
    slopeY, interceptY = np.polyfit(x_ffy, y, 1)
    slopeX, interceptX = np.polyfit(x_ffx, y, 1)

    # Save models to JSON
    model_dict = {
        "FFY": {"slope": float(slopeY), "intercept": float(interceptY)},
        "FFX": {"slope": float(slopeX), "intercept": float(interceptX)}
    }
    with open(model_out, 'w') as outjson:
        json.dump(model_dict, outjson, indent=2)


    # --------------------------------------------------
    # Write evaluation CSV
    # --------------------------------------------------
    out_csv.write(
        "bed_file,fetal_fraction,predicted_fetal_fraction,chrY_median_fraction\n"
    )

    mae_ffx=[]
    mae_ffy=[]
    ffy_list=[]
    ffx_list=[]
    ff_list=[]

    for bed_path, ff in samples:
        median_frac_y = per_sample_median_y[bed_path]
        predicted_ffy = slopeY * median_frac_y + interceptY
        median_frac_x = per_sample_median_x[bed_path]
        predicted_ffx = slopeX * median_frac_x + interceptX

        ffy_list.append(predicted_ffy)
        ffx_list.append(predicted_ffx)
        ff_list.append(ff)

        mae_ffy.append(abs(predicted_ffy-ff))
        mae_ffx.append(abs(predicted_ffx-ff))

        out_csv.write(
            f"{bed_path},{ff:.6f},{predicted_ffy:.6f},{predicted_ffx:.6f},{median_frac_y:.6f},{median_frac_x:.6f}\n"
        )


    out_csv.write("\n")
    corr_matrix = np.corrcoef(ff_list,ffy_list)
    R_sq = corr_matrix[0, 1] ** 2
    out_csv.write(f"Mean absolute error FFY:{np.average(mae_ffy)} , R2:{R_sq}\n")

    corr_matrix = np.corrcoef(ff_list,ffx_list)
    R_sq = corr_matrix[0, 1] ** 2
    out_csv.write(f"Mean absolute error FFX:{np.average(mae_ffx)} , R2:{R_sq}\n")

