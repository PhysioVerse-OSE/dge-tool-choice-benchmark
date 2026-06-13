from pathlib import Path
import argparse
import re
import pandas as pd
import numpy as np


def group_cols(columns, group):
    pattern = re.compile(rf"^{re.escape(group)} \([0-9]+\)$")
    return [c for c in columns if pattern.match(c)]


def write_subset(df, cols, out_path):
    out = df[[df.columns[0]] + cols]
    out.to_csv(out_path, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="GSE196134-Ready-45.csv")
    parser.add_argument("--outdir", default="figure1_inputs")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--subsample-reps", type=int, default=20)
    parser.add_argument("--swap-reps", type=int, default=20)
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    subsample_dir = outdir / "subsamples"
    swap_dir = outdir / "swaps"
    subsample_dir.mkdir(parents=True, exist_ok=True)
    swap_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path)
    gene_col = df.columns[0]
    if gene_col != "Gene Symbol":
        df = df.rename(columns={gene_col: "Gene Symbol"})

    control_cols = group_cols(df.columns, "Control")
    rsvb_cols = group_cols(df.columns, "RSVB")
    if len(control_cols) != 45 or len(rsvb_cols) != 45:
        raise ValueError(f"Expected 45 Control and 45 RSVB samples, found {len(control_cols)} and {len(rsvb_cols)}.")

    rng = np.random.default_rng(args.seed)

    write_subset(df, control_cols + rsvb_cols, subsample_dir / "subsample_n45_rep01.csv")
    write_subset(df, control_cols + rsvb_cols, swap_dir / "swap0_rep01.csv")

    manifest_rows = []
    manifest_rows.append({"analysis": "subsample", "n_per_group": 45, "replicate": 1,
                          "swap_num": 0, "file": str(subsample_dir / "subsample_n45_rep01.csv")})

    for n in [5, 10, 20]:
        for rep in range(1, args.subsample_reps + 1):
            selected_control = rng.choice(control_cols, size=n, replace=False).tolist()
            selected_rsvb = rng.choice(rsvb_cols, size=n, replace=False).tolist()
            out_file = subsample_dir / f"subsample_n{n}_rep{rep:02d}.csv"
            write_subset(df, selected_control + selected_rsvb, out_file)
            manifest_rows.append({"analysis": "subsample", "n_per_group": n, "replicate": rep,
                                  "swap_num": 0, "file": str(out_file)})

    manifest_rows.append({"analysis": "swap", "n_per_group": 45, "replicate": 1,
                          "swap_num": 0, "file": str(swap_dir / "swap0_rep01.csv")})

    for swap_num in [1, 2, 3, 4, 5]:
        for rep in range(1, args.swap_reps + 1):
            swapped = df.copy()
            selected_control = rng.choice(control_cols, size=swap_num, replace=False).tolist()
            selected_rsvb = rng.choice(rsvb_cols, size=swap_num, replace=False).tolist()
            for c_col, r_col in zip(selected_control, selected_rsvb):
                c_values = swapped[c_col].copy()
                swapped[c_col] = swapped[r_col]
                swapped[r_col] = c_values
            out_file = swap_dir / f"swap{swap_num}_rep{rep:02d}.csv"
            swapped.to_csv(out_file, index=False)
            manifest_rows.append({"analysis": "swap", "n_per_group": 45, "replicate": rep,
                                  "swap_num": swap_num, "file": str(out_file)})

    pd.DataFrame(manifest_rows).to_csv(outdir / "manifest.csv", index=False)
    print(f"Wrote inputs and manifest to {outdir}")


if __name__ == "__main__":
    main()
