#!/usr/bin/env python3
"""
Create one balanced combined figure for SARS-CoV-2 cross-study gene selection.

Panel (a):
    Compact grouped bar plot showing reproducible training genes for edgeR and DESeq2.

Panels (b) to (e):
    Venn-style panels showing edgeR-unique, shared, and DESeq2-unique reproducible genes
    for each leave-one-dataset-out fold.

Expected input result files:
    EdgeR-GSE152418.csv
    EdgeR-GSE161731.csv
    EdgeR-GSE171110.csv
    EdgeR-PMC8202013.csv
    DESeq2-GSE152418.csv
    DESeq2-GSE161731.csv
    DESeq2-GSE171110.csv
    DESeq2-PMC8202013.csv

Each file must contain:
    GeneSymbol
    logFC
    FDR
"""

from pathlib import Path
import argparse
import math

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Circle
import numpy as np


DATASETS = ["GSE152418", "GSE161731", "GSE171110", "PMC8202013"]


# Global font sizes for consistent figure styling
PANEL_LABEL_SIZE = 13
TITLE_SIZE = 11.5
AXIS_LABEL_SIZE = 12
TICK_LABEL_SIZE = 10
LEGEND_SIZE = 10.5
COUNT_LABEL_SIZE = 13
VENN_REGION_SIZE = 11.5
VENN_TOTAL_SIZE = 10.5


def find_result_file(results_dir, tool, dataset):
    candidates = [
        results_dir / f"{tool}-{dataset}.csv",
        results_dir / f"{tool}_{dataset}.csv",
        results_dir / f"{tool}_{dataset}_results.csv",
        results_dir / f"{tool}-{dataset}_results.csv",
    ]

    for path in candidates:
        if path.exists():
            return path

    matches = list(results_dir.glob(f"*{tool}*{dataset}*.csv"))
    if len(matches) == 1:
        return matches[0]

    if len(matches) > 1:
        raise FileNotFoundError(
            f"Multiple files matched {tool} and {dataset}: "
            + ", ".join(str(p) for p in matches)
        )

    raise FileNotFoundError(
        f"No result file found for tool={tool}, dataset={dataset} in {results_dir}"
    )


def load_sig_genes(path, fdr_cut=0.05, logfc_cut=1.0):
    df = pd.read_csv(path)

    required = {"GeneSymbol", "FDR", "logFC"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")

    df = df.dropna(subset=["GeneSymbol", "FDR", "logFC"])
    sig = df[(df["FDR"] < fdr_cut) & (df["logFC"].abs() > logfc_cut)]
    return set(sig["GeneSymbol"].astype(str))


def get_reproducible_training_genes(all_sig, tool, train_sets):
    gene_sets = [all_sig[(tool, dataset)] for dataset in train_sets]
    return set.intersection(*gene_sets)


def format_training_sets(training_text):
    return training_text.replace(";", ", ")


def build_summary(results_dir, fdr_cut=0.05, logfc_cut=1.0):
    all_sig = {}
    input_files = []

    for tool in ["EdgeR", "DESeq2"]:
        for dataset in DATASETS:
            path = find_result_file(results_dir, tool, dataset)
            genes = load_sig_genes(path, fdr_cut=fdr_cut, logfc_cut=logfc_cut)
            all_sig[(tool, dataset)] = genes
            input_files.append(
                {
                    "tool": tool,
                    "dataset": dataset,
                    "file": str(path),
                    "n_sig_genes": len(genes),
                }
            )

    input_df = pd.DataFrame(input_files)

    rows = []
    gene_lists = {}

    for held_out in DATASETS:
        train_sets = [d for d in DATASETS if d != held_out]

        edger_train = get_reproducible_training_genes(all_sig, "EdgeR", train_sets)
        deseq_train = get_reproducible_training_genes(all_sig, "DESeq2", train_sets)

        shared = edger_train & deseq_train
        edger_unique = edger_train - deseq_train
        deseq_unique = deseq_train - edger_train
        union = edger_train | deseq_train

        gene_lists[held_out] = {
            "edgeR_reproducible": edger_train,
            "DESeq2_reproducible": deseq_train,
            "shared": shared,
            "edgeR_unique": edger_unique,
            "DESeq2_unique": deseq_unique,
        }

        rows.append(
            {
                "held_out_dataset": held_out,
                "training_datasets": ";".join(train_sets),
                "edgeR_reproducible_training": len(edger_train),
                "DESeq2_reproducible_training": len(deseq_train),
                "shared_reproducible": len(shared),
                "edgeR_unique": len(edger_unique),
                "DESeq2_unique": len(deseq_unique),
                "union_reproducible": len(union),
                "jaccard_edgeR_DESeq2": len(shared) / len(union) if len(union) > 0 else math.nan,
                "DO_edgeR_vs_DESeq2": len(shared) / len(deseq_train) if len(deseq_train) > 0 else math.nan,
                "DO_DESeq2_vs_edgeR": len(shared) / len(edger_train) if len(edger_train) > 0 else math.nan,
            }
        )

    summary_df = pd.DataFrame(rows)
    return input_df, summary_df, gene_lists


def save_gene_lists(gene_lists, outdir):
    for held_out, lists in gene_lists.items():
        pd.DataFrame({"GeneSymbol": sorted(lists["edgeR_reproducible"])}).to_csv(
            outdir / f"{held_out}_edgeR_reproducible_training_genes.csv", index=False
        )
        pd.DataFrame({"GeneSymbol": sorted(lists["DESeq2_reproducible"])}).to_csv(
            outdir / f"{held_out}_DESeq2_reproducible_training_genes.csv", index=False
        )
        pd.DataFrame({"GeneSymbol": sorted(lists["shared"])}).to_csv(
            outdir / f"{held_out}_shared_reproducible_genes.csv", index=False
        )
        pd.DataFrame({"GeneSymbol": sorted(lists["edgeR_unique"])}).to_csv(
            outdir / f"{held_out}_edgeR_unique_training_signature.csv", index=False
        )
        pd.DataFrame({"GeneSymbol": sorted(lists["DESeq2_unique"])}).to_csv(
            outdir / f"{held_out}_DESeq2_unique_training_signature.csv", index=False
        )


def draw_panel_a_counts(ax, summary_df):
    x = np.arange(len(summary_df))
    width = 0.34

    edge_vals = summary_df["edgeR_reproducible_training"].to_numpy()
    deseq_vals = summary_df["DESeq2_reproducible_training"].to_numpy()

    edge_color = "#2C7FB8"
    deseq_color = "#E76F51"

    bars1 = ax.bar(
        x - width / 2,
        edge_vals,
        width,
        label="edgeR",
        color=edge_color,
        edgecolor="black",
        linewidth=1.1,
    )

    bars2 = ax.bar(
        x + width / 2,
        deseq_vals,
        width,
        label="DESeq2",
        color=deseq_color,
        edgecolor="black",
        linewidth=1.1,
    )

    ymax = max(max(edge_vals), max(deseq_vals))
    ax.set_ylim(0, ymax * 1.23)

    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + ymax * 0.025,
                f"{int(height)}",
                ha="center",
                va="bottom",
                fontsize=COUNT_LABEL_SIZE,
                fontweight="normal",
            )

    x_labels = []
    for _, row in summary_df.iterrows():
        held = row["held_out_dataset"]
        x_labels.append(f"Held-out\n{held}")

    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=TICK_LABEL_SIZE, fontweight="normal")
    ax.set_ylabel("Reproducible\ntraining genes", fontsize=AXIS_LABEL_SIZE, fontweight="normal")
    ax.set_xlabel("Leave-one-dataset-out fold", fontsize=AXIS_LABEL_SIZE, fontweight="normal")

    ax.legend(
        frameon=True,
        fontsize=LEGEND_SIZE,
        loc="upper right",
        ncol=2,
        edgecolor="black",
        fancybox=False,
        framealpha=1.0,
    )

    ax.grid(axis="y", alpha=0.25)

    ax.tick_params(axis="y", labelsize=TICK_LABEL_SIZE)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    ax.text(
        -0.045, 1.035, "(a)", transform=ax.transAxes,
        fontsize=PANEL_LABEL_SIZE, fontweight="bold", va="top", ha="right"
    )


def draw_custom_venn_panel(ax, row, panel_label):
    edge_total = int(row["edgeR_reproducible_training"])
    deseq_total = int(row["DESeq2_reproducible_training"])
    shared = int(row["shared_reproducible"])
    edge_unique = int(row["edgeR_unique"])
    deseq_unique = int(row["DESeq2_unique"])

    held_out = row["held_out_dataset"]
    training_sets = format_training_sets(row["training_datasets"])

    edge_color = "#2C7FB8"
    deseq_color = "#E76F51"
    shared_color = "#8E8E8E"

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left_circle = Circle((0.43, 0.55), 0.25, color=edge_color, alpha=0.45, ec="black", lw=1.2)
    right_circle = Circle((0.57, 0.55), 0.25, color=deseq_color, alpha=0.45, ec="black", lw=1.2)

    ax.add_patch(left_circle)
    ax.add_patch(right_circle)

    overlap_circle = Circle((0.50, 0.55), 0.10, color=shared_color, alpha=0.30, ec=None)
    ax.add_patch(overlap_circle)

    ax.text(
        0.50,
        0.55,
        f"{shared}\nshared",
        ha="center",
        va="center",
        fontsize=VENN_REGION_SIZE,
        fontweight="normal",
        color="black",
    )

    ax.annotate(
        f"{edge_unique}\nedgeR-unique",
        xy=(0.25, 0.55),
        xytext=(0.09, 0.55),
        ha="center",
        va="center",
        fontsize=VENN_REGION_SIZE,
        fontweight="normal",
        arrowprops=dict(arrowstyle="-", lw=1.1, color="black"),
    )

    ax.annotate(
        f"{deseq_unique}\nDESeq2-unique",
        xy=(0.75, 0.55),
        xytext=(0.91, 0.55),
        ha="center",
        va="center",
        fontsize=VENN_REGION_SIZE,
        fontweight="normal",
        arrowprops=dict(arrowstyle="-", lw=1.1, color="black"),
    )

    ax.text(
        0.30,
        0.82,
        f"edgeR reproducible\nn={edge_total}",
        ha="center",
        va="center",
        fontsize=VENN_TOTAL_SIZE,
        fontweight="normal",
        color="black",
    )

    ax.text(
        0.70,
        0.82,
        f"DESeq2 reproducible\nn={deseq_total}",
        ha="center",
        va="center",
        fontsize=VENN_TOTAL_SIZE,
        fontweight="normal",
        color="black",
    )

    ax.text(
        0.02,
        0.98,
        panel_label,
        transform=ax.transAxes,
        fontsize=PANEL_LABEL_SIZE,
        fontweight="bold",
        va="top",
        ha="left",
    )

    ax.text(
        0.50,
        0.08,
        f"Held-out: {held_out}\nTraining: {training_sets}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=TITLE_SIZE,
        fontweight="normal",
    )


def make_single_combined_figure(summary_df, outdir):
    fig = plt.figure(figsize=(12, 12.2))

    gs = fig.add_gridspec(
        4, 2,
        height_ratios=[0.52, 0.13, 1.0, 1.0],
        hspace=0.18,
        wspace=0.14
    )

    ax_a = fig.add_subplot(gs[0, :])
    draw_panel_a_counts(ax_a, summary_df)

    legend_ax = fig.add_subplot(gs[1, :])
    legend_ax.axis("off")

    handles = [
        Patch(facecolor="#2C7FB8", edgecolor="black", alpha=0.45, label="edgeR reproducible training genes"),
        Patch(facecolor="#E76F51", edgecolor="black", alpha=0.45, label="DESeq2 reproducible training genes"),
        Patch(facecolor="#8E8E8E", edgecolor="black", alpha=0.30, label="Shared between tools"),
    ]

    legend = legend_ax.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=3,
        frameon=True,
        fontsize=LEGEND_SIZE,
        fancybox=False,
        edgecolor="black",
        framealpha=1.0,
        handlelength=1.5,
        columnspacing=1.8,
        borderpad=0.6,
    )

    legend.get_frame().set_linewidth(1.0)

    venn_axes = [
        fig.add_subplot(gs[2, 0]),
        fig.add_subplot(gs[2, 1]),
        fig.add_subplot(gs[3, 0]),
        fig.add_subplot(gs[3, 1]),
    ]

    panel_labels = ["(b)", "(c)", "(d)", "(e)"]

    for ax, (_, row), label in zip(venn_axes, summary_df.iterrows(), panel_labels):
        draw_custom_venn_panel(ax, row, label)

    fig.tight_layout(rect=[0.025, 0.025, 0.985, 0.99])

    fig.savefig(outdir / "Figure_S1_combined_training_gene_selection_balanced.png", dpi=600, bbox_inches="tight")
    fig.savefig(outdir / "Figure_S1_combined_training_gene_selection_balanced.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-dir",
        default=".",
        help="Folder containing edgeR and DESeq2 result CSV files.",
    )
    parser.add_argument(
        "--outdir",
        default="cross_study_venn_outputs",
        help="Output folder.",
    )
    parser.add_argument("--fdr-cut", type=float, default=0.05)
    parser.add_argument("--logfc-cut", type=float, default=1.0)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_df, summary_df, gene_lists = build_summary(
        results_dir=results_dir,
        fdr_cut=args.fdr_cut,
        logfc_cut=args.logfc_cut,
    )

    input_df.to_csv(outdir / "input_deg_file_summary.csv", index=False)
    summary_df.to_csv(outdir / "cross_study_training_venn_summary.csv", index=False)
    save_gene_lists(gene_lists, outdir)

    make_single_combined_figure(summary_df, outdir)

    print(f"Wrote outputs to {outdir}")
    print(summary_df)


if __name__ == "__main__":
    main()