from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


def summarize_directional_overlap(do_df, n_order):
    rows = []
    for n in n_order:
        sub = do_df[do_df["n_per_group"] == n]
        for col in ["DO_edgeR_vs_DESeq2", "DO_DESeq2_vs_edgeR"]:
            vals = sub[col].dropna().to_numpy(dtype=float)
            mean_val = np.mean(vals) if len(vals) > 0 else np.nan
            sd_val = np.std(vals, ddof=1) if len(vals) > 1 else 0.0
            rows.append(
                {
                    "n_per_group": n,
                    "metric": col,
                    "mean": mean_val,
                    "sd": sd_val,
                    "n_reps": len(vals),
                }
            )
    return pd.DataFrame(rows)


def draw_panel_i(ax, do_df):
    n_order = [5, 10, 20, 45]
    summary = summarize_directional_overlap(do_df, n_order)

    purple = "#800080"
    green = "#006400"

    x = np.array(n_order, dtype=float)

    y1 = summary.loc[
        summary["metric"] == "DO_edgeR_vs_DESeq2", "mean"
    ].to_numpy(dtype=float)
    e1 = summary.loc[
        summary["metric"] == "DO_edgeR_vs_DESeq2", "sd"
    ].to_numpy(dtype=float)

    y2 = summary.loc[
        summary["metric"] == "DO_DESeq2_vs_edgeR", "mean"
    ].to_numpy(dtype=float)
    e2 = summary.loc[
        summary["metric"] == "DO_DESeq2_vs_edgeR", "sd"
    ].to_numpy(dtype=float)

    print("Panel i purple means:", y1)
    print("Panel i purple SDs:", e1)
    print("Panel i green means:", y2)
    print("Panel i green SDs:", e2)

    for xv in x:
        ax.axvline(x=xv, color="#bfbfbf", lw=2, alpha=0.8, zorder=0)

    ax.grid(axis="y", linestyle="--", color="#bfbfbf", alpha=0.8)

    ax.errorbar(
        x, y1, yerr=e1,
        fmt="o-",
        color=purple,
        lw=3.2,
        ms=14,
        capsize=8,
        capthick=2.2,
        elinewidth=2.2,
        markeredgecolor="white",
        markeredgewidth=1.2,
        label="DO(edgeR, DESeq2)",
        zorder=3
    )

    ax.errorbar(
        x, y2, yerr=e2,
        fmt="s-",
        color=green,
        lw=3.2,
        ms=13,
        capsize=8,
        capthick=2.2,
        elinewidth=2.2,
        markeredgecolor="white",
        markeredgewidth=1.2,
        label="DO(DESeq2, edgeR)",
        zorder=4
    )

    for xv, yv, ev in zip(x, y1, e1):
        ax.text(
            xv, yv - 0.10, f"{yv:.2f} ± {ev:.2f}",
            color=purple, fontsize=14, ha="center", va="top"
        )

    for xv, yv, ev in zip(x, y2, e2):
        ax.text(
            xv, yv + 0.04, f"{yv:.2f} ± {ev:.2f}",
            color=green, fontsize=14, ha="center", va="bottom"
        )

    ax.set_xlim(3, 47)
    ax.set_ylim(0.0, 1.05)
    ax.set_xticks(n_order)
    ax.set_xticklabels([str(n) for n in n_order], fontsize=18)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.tick_params(axis="y", labelsize=18)

    ax.set_xlabel("Sample Size", fontsize=22, fontweight="bold")
    ax.set_ylabel("Directional Overlap (DO)", fontsize=22, fontweight="bold")

    leg = ax.legend(
        title="Comparison",
        fontsize=16,
        title_fontsize=18,
        loc="lower right",
        frameon=True
    )
    leg.get_frame().set_edgecolor("#bfbfbf")
    leg.get_frame().set_linewidth(2)

    ax.text(
        -0.08, 0.5, "(i)", transform=ax.transAxes,
        fontsize=20, fontweight="bold", va="center", ha="right"
    )

    for spine in ax.spines.values():
        spine.set_linewidth(2.5)
        spine.set_color("#bfbfbf")

def draw_panel_j(ax, jac_df):
    rng = np.random.default_rng(12345)

    swap_order = [1, 2, 3, 4, 5]
    edge_color = "#1f77b4"
    deseq_color = "#e07b7b"

    pos_edge = np.array(swap_order, dtype=float) - 0.16
    pos_deseq = np.array(swap_order, dtype=float) + 0.16

    edge_data = [
        jac_df[(jac_df["tool"] == "edgeR") & (jac_df["swap_num"] == s)]["jaccard_vs_unperturbed"]
        .dropna()
        .to_numpy(dtype=float)
        for s in swap_order
    ]
    deseq_data = [
        jac_df[(jac_df["tool"] == "DESeq2") & (jac_df["swap_num"] == s)]["jaccard_vs_unperturbed"]
        .dropna()
        .to_numpy(dtype=float)
        for s in swap_order
    ]

    bp1 = ax.boxplot(
        edge_data,
        positions=pos_edge,
        widths=0.28,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.8),
        whiskerprops=dict(color="#666666", linewidth=1.6),
        capprops=dict(color="#666666", linewidth=1.6),
        boxprops=dict(color=edge_color, linewidth=1.6),
    )
    for patch in bp1["boxes"]:
        patch.set_facecolor(edge_color)
        patch.set_alpha(0.85)

    bp2 = ax.boxplot(
        deseq_data,
        positions=pos_deseq,
        widths=0.28,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.8),
        whiskerprops=dict(color="#666666", linewidth=1.6),
        capprops=dict(color="#666666", linewidth=1.6),
        boxprops=dict(color=deseq_color, linewidth=1.6),
    )
    for patch in bp2["boxes"]:
        patch.set_facecolor(deseq_color)
        patch.set_alpha(0.85)

    for p, vals in zip(pos_edge, edge_data):
        jitter = rng.uniform(-0.065, 0.065, size=len(vals))
        ax.scatter(
            np.full(len(vals), p) + jitter,
            vals,
            s=38,
            color="blue",
            edgecolor="black",
            linewidth=0.45,
            alpha=0.8,
            zorder=3,
        )

    for p, vals in zip(pos_deseq, deseq_data):
        jitter = rng.uniform(-0.065, 0.065, size=len(vals))
        ax.scatter(
            np.full(len(vals), p) + jitter,
            vals,
            s=38,
            color="red",
            edgecolor="black",
            linewidth=0.45,
            alpha=0.8,
            zorder=3,
        )

    for xv in [1.5, 2.5, 3.5, 4.5]:
        ax.axvline(x=xv, color="#808080", linestyle="--", linewidth=1.5, alpha=0.8, zorder=0)

    ax.grid(axis="y", linestyle="-", color="#bfbfbf", linewidth=1.5, alpha=0.8)

    ax.set_xlim(0.5, 5.5)
    ax.set_ylim(0.40, 1.00)
    ax.set_xticks(swap_order)
    ax.set_xticklabels([str(s) for s in swap_order], fontsize=18)
    ax.tick_params(axis="y", labelsize=18)

    ax.set_xlabel("Number of Swapped Samples (Outliers)", fontsize=22, fontweight="bold")
    ax.set_ylabel("Jaccard Index", fontsize=22, fontweight="bold")

    legend_elements = [
        Patch(facecolor=edge_color, edgecolor=edge_color, alpha=0.85, label="edgeR"),
        Patch(facecolor=deseq_color, edgecolor=deseq_color, alpha=0.85, label="DESeq2"),
    ]
    leg = ax.legend(
        handles=legend_elements,
        title="Tool",
        fontsize=18,
        title_fontsize=20,
        loc="upper right",
        ncol=2,
        frameon=True,
        borderpad=0.8,
    )
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(1.8)

    ax.text(
        -0.08, 0.5, "(j)", transform=ax.transAxes,
        fontsize=20, fontweight="bold", va="center", ha="right"
    )

    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
        spine.set_color("#bfbfbf")

def save_separate_panel_i(summary_dir, outdir):
    do_df = pd.read_csv(summary_dir / "directional_overlap_by_replicate.csv")
    fig, ax = plt.subplots(figsize=(14, 5))
    draw_panel_i(ax, do_df)
    fig.tight_layout()
    fig.savefig(outdir / "panel_i_directional_overlap.png", dpi=600, bbox_inches="tight")
    fig.savefig(outdir / "panel_i_directional_overlap.pdf", bbox_inches="tight")
    plt.close(fig)


def save_separate_panel_j(summary_dir, outdir):
    jac_df = pd.read_csv(summary_dir / "outlier_jaccard_by_replicate.csv")
    fig, ax = plt.subplots(figsize=(14, 5))
    draw_panel_j(ax, jac_df)
    fig.tight_layout()
    fig.savefig(outdir / "panel_j_outlier_jaccard.png", dpi=600, bbox_inches="tight")
    fig.savefig(outdir / "panel_j_outlier_jaccard.pdf", bbox_inches="tight")
    plt.close(fig)


def save_combined_figure(summary_dir, outdir):
    do_df = pd.read_csv(summary_dir / "directional_overlap_by_replicate.csv")
    jac_df = pd.read_csv(summary_dir / "outlier_jaccard_by_replicate.csv")

    fig, axes = plt.subplots(
        2, 1, figsize=(14, 10),
        gridspec_kw={"height_ratios": [1.0, 1.0]}
    )

    draw_panel_i(axes[0], do_df)
    draw_panel_j(axes[1], jac_df)

    fig.tight_layout(h_pad=2.0)
    fig.savefig(outdir / "Figure1_panels_i_j_revised.png", dpi=600, bbox_inches="tight")
    fig.savefig(outdir / "Figure1_panels_i_j_revised.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-dir", default="figure1_summary")
    parser.add_argument("--outdir", default="figure1_plots")
    args = parser.parse_args()

    summary_dir = Path(args.summary_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    save_separate_panel_i(summary_dir, outdir)
    save_separate_panel_j(summary_dir, outdir)
    save_combined_figure(summary_dir, outdir)

    print(f"Wrote manuscript-style plots to {outdir}")


if __name__ == "__main__":
    main()