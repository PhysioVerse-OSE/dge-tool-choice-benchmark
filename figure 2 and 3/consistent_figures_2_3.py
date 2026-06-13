#!/usr/bin/env python3
"""
Regenerate Figures 2 and 3 with consistent dataset/contrast labels.

All labels use the format:
    GEO/PMC accession ID + biological contrast

This script replaces ambiguous labels such as:
    Influenza data1
    Influenza data2
    Fibrosis
with:
    PMC8202013 | Influenza
    GSE161731 | Influenza
    GSE231693 | IPF
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


OUTDIR = Path("consistent_figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Central label dictionary
# ---------------------------------------------------------------------

keys = [
    "MPXV-DPI3", "MPXV-DPI7", "MPXV-DPI10", "MPXV-DPI14",
    "EBOV-DPI3", "EBOV-DPI5", "EBOV-DPI7", "EBOV-NEC",
    "SARS-CoV-2", "PMC8202013-Influenza",
    "Bacterial", "GSE161731-Influenza",
    "IPF"
]

labels_long = {
    "MPXV-DPI3": "GSE234118 | Mpox DPI 3",
    "MPXV-DPI7": "GSE234118 | Mpox DPI 7",
    "MPXV-DPI10": "GSE234118 | Mpox DPI 10",
    "MPXV-DPI14": "GSE234118 | Mpox DPI 14",
    "EBOV-DPI3": "GSE115785 | EBOV DPI 3",
    "EBOV-DPI5": "GSE115785 | EBOV DPI 5",
    "EBOV-DPI7": "GSE115785 | EBOV DPI 7",
    "EBOV-NEC": "GSE115785 | EBOV NEC",
    "SARS-CoV-2": "PMC8202013 | SARS-CoV-2",
    "PMC8202013-Influenza": "PMC8202013 | Influenza",
    "Bacterial": "GSE161731 | Bacterial",
    "GSE161731-Influenza": "GSE161731 | Influenza",
    "IPF": "GSE231693 | IPF",
}

labels_multiline = {
    k: labels_long[k].replace(" | ", "\n") for k in keys
}

# ---------------------------------------------------------------------
# Figure 2 data
# ---------------------------------------------------------------------

edger_up = np.array([34, 50, 32, 39, 0, 7, 107, 37, 9, 71, 73, 25, 17])
deseq2_up = np.array([19, 37, 28, 13, 36, 129, 908, 1439, 107, 35, 153, 101, 317])
edger_down = np.array([1, 5, 2, 2, 9, 20, 667, 762, 29, 40, 81, 19, 13])
deseq2_down = np.array([56, 253, 77, 10, 0, 70, 28, 12, 5, 59, 246, 44, 164])

up_jaccard = np.array([0.876, 0.842, 0.845, 0.804, 0.163, 0.857, 0.785, 0.698, 0.897, 0.932, 0.862, 0.639, 0.651])
down_jaccard = np.array([0.313, 0.494, 0.470, 0.455, 0.100, 0.408, 0.372, 0.422, 0.832, 0.876, 0.771, 0.640, 0.665])

pearson_corr = np.array([0.354, 0.283, 0.329, 0.504, 0.027, 0.246, 0.093, 0.094, 0.300, 0.384, 0.284, 0.194, 0.275])
spearman_corr = np.array([0.867, 0.813, 0.783, 0.788, 0.690, 0.898, 0.412, 0.535, 0.870, 0.840, 0.862, 0.623, 0.810])


def make_figure2():
    x = np.arange(len(keys))
    xlabels = [labels_multiline[k] for k in keys]

    fig, axes = plt.subplots(3, 1, figsize=(14, 15), sharex=False)

    # Panel a
    ax = axes[0]
    ax.plot(x, np.log2(edger_up + 1), marker="o", linewidth=2.5, label="edgeR Upregulated", color="#1f77b4")
    ax.plot(x, np.log2(deseq2_up + 1), marker="s", linewidth=2.5, label="DESeq2 Upregulated", color="#ff7f0e")
    ax.plot(x, np.log2(edger_down + 1), marker="^", linestyle="--", linewidth=2.5, label="edgeR Downregulated", color="#2ca02c")
    ax.plot(x, np.log2(deseq2_down + 1), marker="v", linestyle="--", linewidth=2.5, label="DESeq2 Downregulated", color="#d62728")
    ax.set_ylabel("log2(Unique Gene Count + 1)", fontsize=15, fontweight="bold")
    ax.legend(fontsize=12, ncol=2, frameon=True)
    ax.grid(True, linestyle="--", alpha=0.45)
    ax.text(0.01, 0.90, "(a)", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")

    # Panel b
    ax = axes[1]
    ax.plot(x, up_jaccard, marker="o", linewidth=2.5, label="Upregulated Genes Jaccard Index", color="#9467bd")
    ax.plot(x, down_jaccard, marker="s", linewidth=2.5, label="Downregulated Genes Jaccard Index", color="#8c564b")
    ax.set_ylabel("Jaccard Index", fontsize=15, fontweight="bold")
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=12, frameon=True)
    ax.grid(True, linestyle="--", alpha=0.45)
    ax.text(0.01, 0.90, "(b)", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")

    # Panel c
    ax = axes[2]
    ax.plot(x, pearson_corr, marker="o", linewidth=2.5, label="Pearson r (Bonferroni p)", color="#e377c2")
    ax.plot(x, spearman_corr, marker="s", linewidth=2.5, label="Spearman ρ (Bonferroni p)", color="#17becf")
    ax.set_ylabel("Correlation Coefficient", fontsize=15, fontweight="bold")
    ax.set_xlabel("Accession ID and biological contrast", fontsize=15, fontweight="bold")
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=12, frameon=True)
    ax.grid(True, linestyle="--", alpha=0.45)
    ax.text(0.01, 0.90, "(c)", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=13)
        ax.tick_params(axis="y", labelsize=12)

    fig.tight_layout(h_pad=1.2)
    fig.savefig(OUTDIR / "Main-2_consistent_labels.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTDIR / "Main-2_consistent_labels.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------
# Figure 3 data
# ---------------------------------------------------------------------

precision_edgeR = np.array([0.882, 1.000, 0.944, 1.000, 0.889, 1.000, 1.000, 1.000, 0.990, 1.000, 1.000, 0.944, 0.952])
precision_DESeq2 = np.array([1.000, 0.867, 1.000, 1.000, 0.833, 0.917, 1.000, 1.000, 0.886, 1.000, 0.913, 0.941, 0.900])
recall_edgeR = np.array([0.833, 0.944, 0.944, 1.000, 0.727, 1.000, 1.000, 1.000, 0.981, 1.000, 1.000, 1.000, 1.000])
recall_DESeq2 = np.array([0.889, 0.722, 1.000, 0.889, 0.909, 1.000, 1.000, 1.000, 0.903, 1.000, 0.913, 0.941, 0.900])
f1_edgeR = np.array([0.857, 0.971, 0.944, 1.000, 0.800, 1.000, 1.000, 1.000, 0.985, 1.000, 1.000, 0.971, 0.976])
f1_DESeq2 = np.array([0.941, 0.788, 1.000, 0.941, 0.870, 0.957, 1.000, 1.000, 0.894, 1.000, 0.913, 0.941, 0.900])


def dolan_more_values(edge_vals, deseq_vals):
    max_vals = np.maximum(edge_vals, deseq_vals)
    ratio_edge = max_vals / edge_vals
    ratio_deseq = max_vals / deseq_vals
    tau = np.linspace(1, max(np.max(ratio_edge), np.max(ratio_deseq)) + 0.01, 200)
    rho_edge = np.array([np.mean(ratio_edge <= t) for t in tau])
    rho_deseq = np.array([np.mean(ratio_deseq <= t) for t in tau])
    return tau, rho_edge, rho_deseq


def radar_axis(ax, edge_vals, deseq_vals, title, panel_label, show_legend=False):
    labels = [labels_multiline[k] for k in keys]
    n = len(labels)

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angles_closed = np.concatenate([angles, angles[:1]])
    edge_closed = np.concatenate([edge_vals, edge_vals[:1]])
    deseq_closed = np.concatenate([deseq_vals, deseq_vals[:1]])

    ax.plot(angles_closed, edge_closed, label="edgeR", linewidth=2.2, marker="o", color="#1f77b4")
    ax.fill(angles_closed, edge_closed, alpha=0.18, color="black")
    ax.plot(angles_closed, deseq_closed, label="DESeq2", linewidth=2.2, marker="s", color="#8B1A1A")
    ax.fill(angles_closed, deseq_closed, alpha=0.15, color="#8B1A1A")

    ax.set_xticks(angles)
    ax.set_xticklabels(labels, fontsize=9.5,fontweight="bold")
    ax.tick_params(axis="x", pad=8)
    ax.set_ylim(0.60, 1.05)
    ax.set_yticks([0.7, 0.8, 0.9, 1.0])
    ax.set_yticklabels(["0.7", "0.8", "0.9", "1.0"], fontsize=8)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=14)
    ax.grid(True, linestyle="--", alpha=0.55)
    ax.text(-0.10, 1.10, panel_label, transform=ax.transAxes, fontsize=16, fontweight="bold")

    if show_legend:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=13, frameon=True)


def make_figure3():
    fig = plt.figure(figsize=(14.5, 14.5))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)

    ax1 = fig.add_subplot(gs[0, 0], polar=True)
    ax2 = fig.add_subplot(gs[0, 1], polar=True)
    ax3 = fig.add_subplot(gs[1, 0], polar=True)
    ax4 = fig.add_subplot(gs[1, 1])

    radar_axis(ax1, precision_edgeR, precision_DESeq2, "Precision Comparison", "(a)")
    radar_axis(ax2, recall_edgeR, recall_DESeq2, "Recall Comparison", "(b)")
    radar_axis(ax3, f1_edgeR, f1_DESeq2, "F1 Score Comparison", "(c)", show_legend=True)

    tau, rho_edge, rho_deseq = dolan_more_values(f1_edgeR, f1_DESeq2)
    ax4.plot(tau, rho_edge, label="edgeR", linewidth=2.6, marker="o", markevery=30, color="blue")
    ax4.plot(tau, rho_deseq, label="DESeq2", linewidth=2.6, marker="s", markevery=30, color="#8B1A1A")
    ax4.set_xlabel(r"$\tau$", fontsize=18, fontweight="bold")
    ax4.set_ylabel(r"$\rho_s(\tau)$", fontsize=18, fontweight="bold")
    ax4.set_title("Dolan-More Performance Profile\n(F1 Score Across All Contrasts)", fontsize=15, fontweight="bold")
    ax4.tick_params(axis="both", labelsize=14)
    ax4.grid(True, linestyle="--", alpha=0.6)
    ax4.legend(fontsize=14, frameon=True)
    ax4.text(-0.10, 1.08, "(d)", transform=ax4.transAxes, fontsize=16, fontweight="bold")

    fig.tight_layout()
    fig.savefig(OUTDIR / "Main-3_consistent_labels.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTDIR / "Main-3_consistent_labels.pdf", bbox_inches="tight")
    plt.close(fig)


def write_label_key():
    rows = ["Original/short key,Consistent label"]
    original_map = {
        "MPXV-DPI3": "MPXV-DPI3",
        "MPXV-DPI7": "MPXV-DPI7",
        "MPXV-DPI10": "MPXV-DPI10",
        "MPXV-DPI14": "MPXV-DPI14",
        "EBOV-DPI3": "EBOV-DPI3",
        "EBOV-DPI5": "EBOV-DPI5",
        "EBOV-DPI7": "EBOV-DPI7",
        "EBOV-DPINEC": "EBOV-NEC",
        "SARS-CoV-2": "SARS-CoV-2",
        "Influenza data1": "PMC8202013-Influenza",
        "Bacterial": "Bacterial",
        "Influenza data2": "GSE161731-Influenza",
        "Fibrosis": "IPF",
    }

    for original, key in original_map.items():
        rows.append(f"{original},{labels_long[key]}")
    (OUTDIR / "consistent_label_key.csv").write_text("\n".join(rows), encoding="utf-8")


def main():
    make_figure2()
    make_figure3()
    write_label_key()
    print(f"Wrote revised figures to: {OUTDIR.resolve()}")


if __name__ == "__main__":
    main()
