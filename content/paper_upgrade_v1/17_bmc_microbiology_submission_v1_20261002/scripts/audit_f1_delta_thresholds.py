#!/usr/bin/env python3
"""Audit the frozen F1 positive-delta count over prespecified thresholds."""
from pathlib import Path
import argparse
import pandas as pd

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=Path(__file__).resolve().parents[1] / "source_data/figure1/07_route_strain_flux_matrix_final.csv")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "source_data/figure1/F1_delta_threshold_sensitivity.tsv")
    args = parser.parse_args()
    frame = pd.read_csv(args.matrix, low_memory=False)
    delta = pd.to_numeric(frame.loc[frame["formal_flux_eligible"].astype(bool), "anaerobic_substrate_attributable_product_delta"], errors="coerce").dropna()
    rows = []
    for threshold in (1e-8, 1e-6, 1e-4):
        rows.append({"threshold_model_units": threshold, "eligible_rows_with_delta": int(len(delta)), "positive_rows": int((delta > threshold).sum()), "min_positive_delta": float(delta.min()), "max_positive_delta": float(delta.max())})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
