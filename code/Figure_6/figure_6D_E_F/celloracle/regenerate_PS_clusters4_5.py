#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pandas as pd

from functions.plotting_utils import plot_ps_clusters4_5


def parse_args():
    parser = argparse.ArgumentParser(description="Regenerate only the HIF1A KO perturbation-score panel.")
    parser.add_argument("--result-dir", default="results/celloracle_hif1a_ko_1500_PS_markov_20260701_133734")
    return parser.parse_args()


def main():
    result_dir = Path(parse_args().result_dir)
    tables, figures = result_dir / "tables", result_dir / "figures"
    grid = pd.read_csv(tables / "HIF1A_KO_perturbation_score_grid.csv", dtype={"seurat_clusters": str})
    cells = pd.read_csv(tables / "HIF1A_KO_markov_density_by_cell.csv", dtype={"seurat_clusters": str})
    output = figures / "HIF1A_KO_perturbation_score_PS_clusters4_5"
    plot_ps_clusters4_5(grid, cells, output)
    print(output.with_suffix(".png").resolve())
    print(output.with_suffix(".pdf").resolve())


if __name__ == "__main__":
    main()
