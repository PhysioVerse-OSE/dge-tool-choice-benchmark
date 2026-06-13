from pathlib import Path
import argparse
import math
import pandas as pd
import numpy as np

try:
    from scipy.stats import wilcoxon
except Exception:
    wilcoxon = None


def sig_genes(path, p_col="Bonferroni_pvalue", logfc_col="logFC", p_cut=0.05, logfc_cut=1.0):
    df = pd.read_csv(path)
    df = df.dropna(subset=[p_col, logfc_col])
    sig = df[(df[p_col] < p_cut) & (df[logfc_col].abs() > logfc_cut)]
    return set(sig["GeneSymbol"].astype(str))


def directional_overlap(a, b):
    if len(b) == 0:
        return math.nan
    return len(a & b) / len(b)


def jaccard(a, b):
    union = a | b
    if len(union) == 0:
        return math.nan
    return len(a & b) / len(union)


def result_path(results_dir, tool, input_file):
    stem = Path(input_file).stem
    prefix = "edgeR" if tool == "edgeR" else "DESeq2"
    return Path(results_dir) / tool / f"{prefix}_{stem}.csv"


def p_value_paired(x, y):
    if wilcoxon is None:
        return np.nan
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or np.allclose(x - y, 0):
        return np.nan
    return wilcoxon(x, y, zero_method="wilcox", alternative="two-sided").pvalue


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="figure1_dge_results")
    parser.add_argument("--outdir", default="figure1_summary")
    parser.add_argument("--p-cut", type=float, default=0.05)
    parser.add_argument("--logfc-cut", type=float, default=1.0)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(results_dir / "manifest.csv")

    gene_sets = {}
    rows = []
    for _, row in manifest.iterrows():
        for tool in ["edgeR", "DESeq2"]:
            p = result_path(results_dir, tool, row["file"])
            genes = sig_genes(p, p_cut=args.p_cut, logfc_cut=args.logfc_cut)
            key = (tool, Path(row["file"]).stem)
            gene_sets[key] = genes
            rows.append({**row.to_dict(), "tool": tool, "input_stem": Path(row["file"]).stem,
                         "deg_count": len(genes)})
    deg_counts = pd.DataFrame(rows)
    deg_counts.to_csv(outdir / "deg_counts_by_replicate.csv", index=False)

    subs = deg_counts[deg_counts["analysis"] == "subsample"].copy()
    count_summary = subs.groupby(["n_per_group", "tool"], as_index=False)["deg_count"].agg(["mean", "std", "median", "min", "max"]).reset_index()
    count_summary.to_csv(outdir / "subsample_deg_count_summary.csv", index=False)

    paired_rows = []
    for n, g in subs[subs["n_per_group"] < 45].groupby("n_per_group"):
        wide = g.pivot(index="replicate", columns="tool", values="deg_count").dropna()
        paired_rows.append({"n_per_group": n,
                            "mean_edgeR": wide["edgeR"].mean(),
                            "mean_DESeq2": wide["DESeq2"].mean(),
                            "mean_difference_DESeq2_minus_edgeR": (wide["DESeq2"] - wide["edgeR"]).mean(),
                            "wilcoxon_p": p_value_paired(wide["DESeq2"], wide["edgeR"]),
                            "n_replicates": len(wide)})
    pd.DataFrame(paired_rows).to_csv(outdir / "paired_tests_deg_counts.csv", index=False)

    do_rows = []
    for _, row in manifest[manifest["analysis"] == "subsample"].iterrows():
        stem = Path(row["file"]).stem
        edger = gene_sets[("edgeR", stem)]
        deseq = gene_sets[("DESeq2", stem)]
        do_rows.append({**row.to_dict(), "input_stem": stem,
                        "DO_edgeR_vs_DESeq2": directional_overlap(edger, deseq),
                        "DO_DESeq2_vs_edgeR": directional_overlap(deseq, edger),
                        "Jaccard_edgeR_DESeq2": jaccard(edger, deseq)})
    do_df = pd.DataFrame(do_rows)
    do_df.to_csv(outdir / "directional_overlap_by_replicate.csv", index=False)
    do_summary = do_df.groupby("n_per_group", as_index=False)[["DO_edgeR_vs_DESeq2", "DO_DESeq2_vs_edgeR", "Jaccard_edgeR_DESeq2"]].agg(["mean", "std", "median", "min", "max"]).reset_index()
    do_summary.to_csv(outdir / "directional_overlap_summary.csv", index=False)

    full_edge = gene_sets[("edgeR", "swap0_rep01")]
    full_deseq = gene_sets[("DESeq2", "swap0_rep01")]
    jac_rows = []
    for _, row in manifest[manifest["analysis"] == "swap"].iterrows():
        if row["swap_num"] == 0:
            continue
        stem = Path(row["file"]).stem
        jac_rows.append({**row.to_dict(), "tool": "edgeR", "input_stem": stem,
                         "jaccard_vs_unperturbed": jaccard(gene_sets[("edgeR", stem)], full_edge)})
        jac_rows.append({**row.to_dict(), "tool": "DESeq2", "input_stem": stem,
                         "jaccard_vs_unperturbed": jaccard(gene_sets[("DESeq2", stem)], full_deseq)})
    jac_df = pd.DataFrame(jac_rows)
    jac_df.to_csv(outdir / "outlier_jaccard_by_replicate.csv", index=False)
    jac_summary = jac_df.groupby(["swap_num", "tool"], as_index=False)["jaccard_vs_unperturbed"].agg(["mean", "std", "median", "min", "max"]).reset_index()
    jac_summary.to_csv(outdir / "outlier_jaccard_summary.csv", index=False)

    swap_tests = []
    for s, g in jac_df.groupby("swap_num"):
        wide = g.pivot(index="replicate", columns="tool", values="jaccard_vs_unperturbed").dropna()
        swap_tests.append({"swap_num": s,
                           "mean_edgeR": wide["edgeR"].mean(),
                           "mean_DESeq2": wide["DESeq2"].mean(),
                           "mean_difference_DESeq2_minus_edgeR": (wide["DESeq2"] - wide["edgeR"]).mean(),
                           "wilcoxon_p": p_value_paired(wide["DESeq2"], wide["edgeR"]),
                           "n_replicates": len(wide)})
    pd.DataFrame(swap_tests).to_csv(outdir / "paired_tests_outlier_jaccard.csv", index=False)
    print(f"Wrote summaries to {outdir}")


if __name__ == "__main__":
    main()
