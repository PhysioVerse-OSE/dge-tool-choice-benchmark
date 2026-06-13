#!/usr/bin/env python3

from pathlib import Path
import argparse
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


COLLECTION_ORDER = ["Hallmark", "KEGG"]

CATEGORY_COLORS = {
    "edgeR-only": "#2C7FB8",
    "Shared": "#7A7A7A",
    "DESeq2-only": "#E76F51",
}


def clean_pathway_name(name, collection):
    text = str(name)
    text = re.sub(r"^HALLMARK_", "", text)
    text = re.sub(r"^KEGG_", "", text)
    text = text.replace("_", " ")
    return text


def clean_contrast_name(contrast):
    return str(contrast).replace("_", " ")


def label_two_lines(dataset, contrast):
    return f"{dataset}\n{clean_contrast_name(contrast)}"


def label_three_lines(dataset, contrast, jaccard):
    if pd.isna(jaccard):
        jtxt = "J=NA"
    else:
        jtxt = f"J={jaccard:.2f}"
    return f"{dataset}\n{clean_contrast_name(contrast)}\n{jtxt}"


def get_collection_rows(summary_df, collection, direction="All"):
    sub = summary_df[
        (summary_df["Collection"] == collection) &
        (summary_df["direction"] == direction)
    ].copy()

    if sub.empty:
        return sub

    sub["label2"] = sub.apply(lambda r: label_two_lines(r["dataset"], r["contrast"]), axis=1)
    sub["label3"] = sub.apply(
        lambda r: label_three_lines(r["dataset"], r["contrast"], r["significant_pathway_jaccard"]),
        axis=1
    )

    return sub


def select_most_divergent(summary_df, collection, max_contrasts=4, direction="All"):
    sub = summary_df[
        (summary_df["Collection"] == collection) &
        (summary_df["direction"] == direction)
    ].copy()

    if sub.empty:
        return pd.DataFrame()

    sub["unique_total"] = (
        sub["edgeR_unique_significant_pathways"].fillna(0) +
        sub["DESeq2_unique_significant_pathways"].fillna(0)
    )
    sub["jaccard_fill"] = sub["significant_pathway_jaccard"].fillna(0)
    sub["union_fill"] = sub["union_significant_pathways"].fillna(0)
    sub["divergence_score"] = sub["unique_total"] + (1 - sub["jaccard_fill"]) * sub["union_fill"]

    sub = sub[sub["union_fill"] > 0].copy()

    sub = sub.sort_values(
        ["divergence_score", "unique_total", "union_fill"],
        ascending=[False, False, False]
    ).head(max_contrasts)

    sub["label2"] = sub.apply(lambda r: label_two_lines(r["dataset"], r["contrast"]), axis=1)
    sub["label3"] = sub.apply(
        lambda r: label_three_lines(r["dataset"], r["contrast"], r["significant_pathway_jaccard"]),
        axis=1
    )

    return sub


def get_significant_by_tool(df, fdr_cut):
    edge = df[(df["tool"] == "edgeR") & (df["FDR"] < fdr_cut)].copy()
    deseq = df[(df["tool"] == "DESeq2") & (df["FDR"] < fdr_cut)].copy()
    return edge, deseq


def extract_selected_pathways(enrichment_df, collection, dataset, contrast, fdr_cut=0.05, top_n=4, direction="All"):
    sub = enrichment_df[
        (enrichment_df["Collection"] == collection) &
        (enrichment_df["dataset"] == dataset) &
        (enrichment_df["contrast"] == contrast) &
        (enrichment_df["direction"] == direction)
    ].copy()

    edge, deseq = get_significant_by_tool(sub, fdr_cut)

    edge_set = set(edge["Pathway"])
    deseq_set = set(deseq["Pathway"])
    shared_set = edge_set & deseq_set
    edge_only_set = edge_set - deseq_set
    deseq_only_set = deseq_set - edge_set

    rows = []

    edge_only = edge[edge["Pathway"].isin(edge_only_set)].copy()
    if not edge_only.empty:
        edge_only = edge_only.sort_values(
            ["FDR", "PValue", "Overlap_count"],
            ascending=[True, True, False]
        ).head(top_n)

        for _, row in edge_only.iterrows():
            rows.append({
                "Collection": collection,
                "dataset": dataset,
                "contrast": contrast,
                "comparison": label_three_lines(dataset, contrast, np.nan),
                "Pathway": row["Pathway"],
                "Pathway_clean": clean_pathway_name(row["Pathway"], collection),
                "Category": "edgeR-only",
                "FDR": row["FDR"],
                "Score": -np.log10(max(row["FDR"], 1e-300))
            })

    shared_edge = edge[edge["Pathway"].isin(shared_set)][["Pathway", "FDR", "PValue", "Overlap_count"]].copy()
    shared_deseq = deseq[deseq["Pathway"].isin(shared_set)][["Pathway", "FDR", "PValue", "Overlap_count"]].copy()

    shared = shared_edge.merge(
        shared_deseq,
        on="Pathway",
        suffixes=("_edgeR", "_DESeq2")
    )

    if not shared.empty:
        shared["Mean_FDR"] = shared[["FDR_edgeR", "FDR_DESeq2"]].mean(axis=1)
        shared["Min_PValue"] = shared[["PValue_edgeR", "PValue_DESeq2"]].min(axis=1)
        shared["Mean_overlap"] = shared[["Overlap_count_edgeR", "Overlap_count_DESeq2"]].mean(axis=1)

        shared = shared.sort_values(
            ["Mean_FDR", "Min_PValue", "Mean_overlap"],
            ascending=[True, True, False]
        ).head(top_n)

        for _, row in shared.iterrows():
            rows.append({
                "Collection": collection,
                "dataset": dataset,
                "contrast": contrast,
                "comparison": label_three_lines(dataset, contrast, np.nan),
                "Pathway": row["Pathway"],
                "Pathway_clean": clean_pathway_name(row["Pathway"], collection),
                "Category": "Shared",
                "FDR": row["Mean_FDR"],
                "Score": -np.log10(max(row["Mean_FDR"], 1e-300))
            })

    deseq_only = deseq[deseq["Pathway"].isin(deseq_only_set)].copy()
    if not deseq_only.empty:
        deseq_only = deseq_only.sort_values(
            ["FDR", "PValue", "Overlap_count"],
            ascending=[True, True, False]
        ).head(top_n)

        for _, row in deseq_only.iterrows():
            rows.append({
                "Collection": collection,
                "dataset": dataset,
                "contrast": contrast,
                "comparison": label_three_lines(dataset, contrast, np.nan),
                "Pathway": row["Pathway"],
                "Pathway_clean": clean_pathway_name(row["Pathway"], collection),
                "Category": "DESeq2-only",
                "FDR": row["FDR"],
                "Score": -np.log10(max(row["FDR"], 1e-300))
            })

    return pd.DataFrame(rows)


def reduce_pathways(plot_df, max_pathways=25):
    if plot_df.empty:
        return plot_df

    keep_blocks = []

    for collection, sub in plot_df.groupby("Collection"):
        rank_df = (
            sub.groupby(["Pathway", "Pathway_clean"])
            .agg(
                n_occurrences=("Pathway", "size"),
                max_score=("Score", "max"),
                min_fdr=("FDR", "min")
            )
            .reset_index()
            .sort_values(["n_occurrences", "max_score", "min_fdr"], ascending=[False, False, True])
        )

        keep = set(rank_df.head(max_pathways)["Pathway"])
        keep_blocks.append(sub[sub["Pathway"].isin(keep)].copy())

    return pd.concat(keep_blocks, ignore_index=True)


def make_figure1_all_counts(summary_df, outdir):
    collections = [c for c in COLLECTION_ORDER if c in set(summary_df["Collection"])]
    if not collections:
        raise ValueError("No collections found in pathway_concordance_summary.csv")

    fig, axes = plt.subplots(len(collections), 1, figsize=(22, 7.5 * len(collections)), squeeze=False)
    axes = axes.flatten()

    for i, collection in enumerate(collections):
        ax = axes[i]
        sub = get_collection_rows(summary_df, collection, direction="All")

        if sub.empty:
            ax.axis("off")
            continue

        x = np.arange(len(sub))
        edge_unique = sub["edgeR_unique_significant_pathways"].to_numpy()
        shared = sub["shared_significant_pathways"].to_numpy()
        deseq_unique = sub["DESeq2_unique_significant_pathways"].to_numpy()
        totals = edge_unique + shared + deseq_unique

        ax.bar(
            x,
            edge_unique,
            color=CATEGORY_COLORS["edgeR-only"],
            edgecolor="black",
            linewidth=1.0,
            label="edgeR-only pathways"
        )
        ax.bar(
            x,
            shared,
            bottom=edge_unique,
            color=CATEGORY_COLORS["Shared"],
            edgecolor="black",
            linewidth=1.0,
            label="Shared pathways"
        )
        ax.bar(
            x,
            deseq_unique,
            bottom=edge_unique + shared,
            color=CATEGORY_COLORS["DESeq2-only"],
            edgecolor="black",
            linewidth=1.0,
            label="DESeq2-only pathways"
        )

        ymax = max(totals) if len(totals) > 0 else 1
        ax.set_ylim(0, ymax * 1.22)

        for j, row in sub.reset_index(drop=True).iterrows():
            if pd.isna(row["significant_pathway_jaccard"]):
                jtxt = "J=NA"
            else:
                jtxt = f"J={row['significant_pathway_jaccard']:.2f}"
            ax.text(
                j,
                totals[j] + ymax * 0.03,
                jtxt,
                ha="center",
                va="bottom",
                fontsize=18,
                fontweight="bold"
            )

        ax.set_xticks(x)
        ax.set_xticklabels(sub["label2"], rotation=0, ha="center", fontsize=16)
        ax.set_ylabel("Number of enriched pathways", fontsize=20, fontweight="bold")
        ax.set_title(collection, fontsize=20, fontweight="bold", pad=12)
        ax.grid(axis="y", alpha=0.25)

        panel_label = "(a)" if i == 0 else "(b)"
        ax.text(-0.06, 1.03, panel_label, transform=ax.transAxes, fontsize=20, fontweight="bold", va="bottom")

    handles = [
        Patch(facecolor=CATEGORY_COLORS["edgeR-only"], edgecolor="black", label="edgeR-only pathways"),
        Patch(facecolor=CATEGORY_COLORS["Shared"], edgecolor="black", label="Shared pathways"),
        Patch(facecolor=CATEGORY_COLORS["DESeq2-only"], edgecolor="black", label="DESeq2-only pathways"),
    ]

    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=True, fontsize=22)
    fig.tight_layout(rect=[0, 0.05, 1, 1])

    fig.savefig(outdir / "Figure1_all_datasets_pathway_counts.png", dpi=600, bbox_inches="tight")
    fig.savefig(outdir / "Figure1_all_datasets_pathway_counts.pdf", bbox_inches="tight")
    plt.close(fig)


def collect_bubble_data(enrichment_df, summary_df, fdr_cut=0.05, top_n=6, max_contrasts=6, max_pathways=15):
    meta_blocks = []
    plot_blocks = []

    for collection in COLLECTION_ORDER:
        selected = select_most_divergent(summary_df, collection, max_contrasts=max_contrasts, direction="All")
        if selected.empty:
            continue

        meta_blocks.append(selected)

        for _, row in selected.iterrows():
            block = extract_selected_pathways(
                enrichment_df=enrichment_df,
                collection=collection,
                dataset=row["dataset"],
                contrast=row["contrast"],
                fdr_cut=fdr_cut,
                top_n=top_n,
                direction="All"
            )
            if not block.empty:
                block["comparison_label"] = label_three_lines(
                    row["dataset"],
                    row["contrast"],
                    row["significant_pathway_jaccard"]
                )
                plot_blocks.append(block)

    if not meta_blocks or not plot_blocks:
        return pd.DataFrame(), pd.DataFrame()

    meta_df = pd.concat(meta_blocks, ignore_index=True)
    plot_df = pd.concat(plot_blocks, ignore_index=True)
    plot_df = reduce_pathways(plot_df, max_pathways=max_pathways)

    return meta_df, plot_df


def make_figure2_selected_bubble(meta_df, plot_df, outdir):
    if meta_df.empty or plot_df.empty:
        raise ValueError("No selected bubble plot data were produced.")

    collections = [c for c in COLLECTION_ORDER if c in set(plot_df["Collection"])]
    if not collections:
        raise ValueError("No collections available for Figure 2.")

    heights = []
    for collection in collections:
        n_pathways = plot_df[plot_df["Collection"] == collection]["Pathway"].nunique()
        heights.append(max(6.5, min(12.5, 0.38 * n_pathways + 2.5)))

    fig, axes = plt.subplots(
        len(collections),
        1,
        figsize=(12, sum(heights)),
        gridspec_kw={"height_ratios": heights},
        squeeze=False
    )
    axes = axes.flatten()

    for i, collection in enumerate(collections):
        ax = axes[i]

        meta_sub = meta_df[meta_df["Collection"] == collection].copy()
        plot_sub = plot_df[plot_df["Collection"] == collection].copy()

        if meta_sub.empty or plot_sub.empty:
            ax.axis("off")
            continue

        comparison_order = meta_sub.apply(
            lambda r: label_three_lines(r["dataset"], r["contrast"], r["significant_pathway_jaccard"]),
            axis=1
        ).tolist()

        pathway_rank = (
            plot_sub.groupby(["Pathway", "Pathway_clean"])
            .agg(
                n_occurrences=("Pathway", "size"),
                max_score=("Score", "max"),
                min_fdr=("FDR", "min")
            )
            .reset_index()
            .sort_values(["n_occurrences", "max_score", "min_fdr"], ascending=[True, True, False])
        )

        pathway_order = pathway_rank["Pathway"].tolist()
        pathway_labels = dict(zip(pathway_rank["Pathway"], pathway_rank["Pathway_clean"]))

        x_map = {label: idx for idx, label in enumerate(comparison_order)}
        y_map = {pathway: idx for idx, pathway in enumerate(pathway_order)}

        score_min = plot_sub["Score"].min()
        score_max = plot_sub["Score"].max()

        def size_scale(series):
            s = series.clip(lower=0, upper=12)
            if score_max == score_min:
                return np.full(len(s), 220.0)
            return 110 + (s - score_min) / (score_max - score_min) * 650

        for category in ["edgeR-only", "Shared", "DESeq2-only"]:
            cat = plot_sub[plot_sub["Category"] == category].copy()
            if cat.empty:
                continue

            ax.scatter(
                [x_map[x] for x in cat["comparison_label"]],
                [y_map[y] for y in cat["Pathway"]],
                s=size_scale(cat["Score"]),
                color=CATEGORY_COLORS[category],
                edgecolor="black",
                linewidth=0.7,
                alpha=0.85
            )

        ax.set_xticks(range(len(comparison_order)))
        ax.set_xticklabels(comparison_order, rotation=0, ha="center", fontsize=16)
        ax.set_yticks(range(len(pathway_order)))
        ax.set_yticklabels([pathway_labels[p] for p in pathway_order], fontsize=12)

        ax.set_title(collection, fontsize=20, fontweight="bold", pad=10)
        ax.set_xlabel("Selected dataset/contrast with largest pathway-level differences", fontsize=15, fontweight="bold")
        ax.set_ylabel("Enriched pathway", fontsize=18, fontweight="bold")
        ax.grid(alpha=0.20)

        panel_label = "(a)" if i == 0 else "(b)"
        ax.text(-0.06, 1.03, panel_label, transform=ax.transAxes, fontsize=20, fontweight="bold", va="bottom")

    category_handles = [
        Line2D([0], [0], marker="o", color="w", label="edgeR-only", markerfacecolor=CATEGORY_COLORS["edgeR-only"], markeredgecolor="black", markersize=10),
        Line2D([0], [0], marker="o", color="w", label="Shared", markerfacecolor=CATEGORY_COLORS["Shared"], markeredgecolor="black", markersize=10),
        Line2D([0], [0], marker="o", color="w", label="DESeq2-only", markerfacecolor=CATEGORY_COLORS["DESeq2-only"], markeredgecolor="black", markersize=10),
    ]

    size_handles = [
        plt.scatter([], [], s=120, color="white", edgecolor="black", label="-log10(FDR) low"),
        plt.scatter([], [], s=350, color="white", edgecolor="black", label="-log10(FDR) medium"),
        plt.scatter([], [], s=700, color="white", edgecolor="black", label="-log10(FDR) high"),
    ]

    fig.legend(
        handles=category_handles + size_handles,
        loc="lower center",
        ncol=3,
        frameon=True,
        fontsize=18
    )

    fig.tight_layout(rect=[0, 0.05, 1, 1])

    fig.savefig(outdir / "Figure2_selected_pathway_difference_bubbleplot.png", dpi=600, bbox_inches="tight")
    fig.savefig(outdir / "Figure2_selected_pathway_difference_bubbleplot.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="pathway_enrichment_outputs_v3")
    parser.add_argument("--outdir", default="final_pathway_figures")
    parser.add_argument("--fdr-cut", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=4)
    parser.add_argument("--max-contrasts", type=int, default=4)
    parser.add_argument("--max-pathways", type=int, default=25)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    enrichment_file = results_dir / "pathway_enrichment_all_results.csv"
    summary_file = results_dir / "pathway_concordance_summary.csv"

    if not enrichment_file.exists():
        raise FileNotFoundError(f"Missing file: {enrichment_file}")
    if not summary_file.exists():
        raise FileNotFoundError(f"Missing file: {summary_file}")

    enrichment_df = pd.read_csv(enrichment_file)
    summary_df = pd.read_csv(summary_file)

    make_figure1_all_counts(summary_df, outdir)

    meta_df, plot_df = collect_bubble_data(
        enrichment_df=enrichment_df,
        summary_df=summary_df,
        fdr_cut=args.fdr_cut,
        top_n=args.top_n,
        max_contrasts=args.max_contrasts,
        max_pathways=args.max_pathways
    )

    meta_df.to_csv(outdir / "Figure2_selected_contrasts.csv", index=False)
    plot_df.to_csv(outdir / "Figure2_selected_pathways.csv", index=False)

    make_figure2_selected_bubble(meta_df, plot_df, outdir)

    print(f"Finished. Output folder: {outdir}")


if __name__ == "__main__":
    main()