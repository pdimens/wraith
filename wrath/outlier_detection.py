#!/usr/bin/env python3
"""
Detect outliers in a barcode-sharing matrix via z-scores + nls prediction bands.
Usage: python outlier_detection.py <input_matrix.csv> <output_prefix> <zscore_threshold> <prediction_level>
Example: python outlier_detection.py matrix.csv results/outlier_detection 3 0.95
"""
import sys
import numpy as np
import polars as pl
from scipy.optimize import curve_fit
from scipy.stats import t as t_dist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

def model(x, a, b, c):
    with np.errstate(over="ignore"):
        return np.exp(np.clip(a + b * np.exp(-x * c), -700, 700))  # ponytail: clip avoids overflow warnings during curve_fit's exploratory steps

def main(argv):
    infile, out_prefix, z_thresh, pred_level = argv[1], argv[2], float(argv[3]), float(argv[4])

    m = np.loadtxt(infile, delimiter=",")
    nrow_idx, ncol_idx = np.triu_indices(m.shape[0], k=1)  # strictly upper triangle, matches R's exclusion of diag+lower
    y = m[nrow_idx, ncol_idx]
    # R used 1-based row/col and 0-based-ish "index" = ncol - nrow (both 1-based there); +1 offset keeps values identical
    nrow = nrow_idx + 1
    ncol = ncol_idx + 1
    x = ncol - nrow

    df = pl.DataFrame({"x": x, "y": y, "nrow": nrow, "ncol": ncol})

    # z-score within each diagonal-distance group, vectorized via polars window expr (ddof=1 to match R's scale())
    df = df.with_columns(((pl.col("y") - pl.col("y").mean().over("x")) / pl.col("y").std(ddof=1).over("x")).alias("z_score"))

    # nls fit: y ~ exp(a + b*exp(-x*c))
    x_np, y_np = df["x"].to_numpy(), df["y"].to_numpy()
    popt, pcov = curve_fit(model, x_np, y_np, p0=[1, 1, 1], maxfev=10000)
    a, b, c = popt
    n, p = len(df), len(popt)
    dof = n - p

    fitted = model(x_np, *popt)
    resid = y_np - fitted
    sigma2 = np.sum(resid ** 2) / dof

    # delta-method gradient of model wrt (a,b,c) at each x, vectorized
    xv = df["x"].to_numpy()
    inner = np.exp(-xv * c)
    f = np.exp(a + b * inner)
    d_a = f
    d_b = f * inner
    d_c = f * b * (-xv) * inner
    J = np.column_stack([d_a, d_b, d_c])  # n x 3

    se_fit2 = np.einsum("ij,jk,ik->i", J, pcov, J)  # diag(J @ pcov @ J.T), vectorized
    se_pred = np.sqrt(se_fit2 + sigma2)

    t_val = t_dist.ppf(1 - (1 - pred_level) / 2, dof)
    q_lo_name, q_hi_name = f"Q{(1 - pred_level) / 2 * 100:.1f}", f"Q{(1 - (1 - pred_level) / 2) * 100:.1f}"
    df = df.with_columns(
        pl.Series("Estimate", fitted),
        pl.Series("Est.Error", np.sqrt(se_fit2)),
        pl.Series("Qbottom", fitted - t_val * se_pred),
        pl.Series("Qtop", fitted + t_val * se_pred),
        pl.col("z_score").abs().alias("abs_z"),
    ).with_columns(
        (((pl.col("y") > pl.col("Qtop")) | (pl.col("y") < pl.col("Qbottom"))) & (pl.col("abs_z") > z_thresh)).alias("is_outlier"),
        (pl.col("y") > pl.col("Qtop")).alias("upper"),
        (pl.col("y") < pl.col("Qbottom")).alias("lower"),
    )

    # plot
    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(6, 3.5))
    sns.scatterplot(x=df["x"], y=df["y"], s=8, color="black", ax=ax)
    order = np.argsort(xv)
    ax.plot(xv[order], fitted[order], color="blue")
    ax.fill_between(xv[order], df["Qbottom"].to_numpy()[order], df["Qtop"].to_numpy()[order],
                     color="purple", alpha=0.3)
    out_df = df.filter(pl.col("is_outlier"))
    sns.scatterplot(x=out_df["x"], y=out_df["y"], color="#FF6600", ax=ax)
    ax.set_xlabel("Distance from matrix")
    ax.set_ylabel("Similarity index")
    ax.set_title(f"{pred_level * 100:g}% prediction bands")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_plot.png", dpi=300)

    outliers = (
        out_df.select(["nrow", "ncol", "y", "Estimate", "Est.Error", "Qbottom", "Qtop", "upper", "lower", "z_score"])
        .rename({"y": "value", "Qbottom": q_lo_name, "Qtop": q_hi_name})
    )
    outliers.write_csv(f"{out_prefix}.csv")


if __name__ == "__main__":
    main(sys.argv)