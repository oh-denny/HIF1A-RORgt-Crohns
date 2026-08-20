#!/usr/bin/env python3

import argparse
from pathlib import Path

from functions.io_utils import configure_runtime_environment

PROJECT_DIR = Path(__file__).resolve().parents[4]
configure_runtime_environment(PROJECT_DIR, include_celloracle_cache=True)

import matplotlib

matplotlib.use("Agg")

from functions.celloracle_utils import load_oracle
from functions.perturbation_utils import celloracle_impact_summary
from functions.plotting_utils import plot_propagation_magnitude


def parse_args():
    parser = argparse.ArgumentParser(
        description="Regenerate the CellOracle HIF1A KO propagation-impact figure."
    )
    parser.add_argument(
        "--oracle",
        default="results/celloracle_hif1a_ko_1500/objects/oracle_HIF1A_KO.celloracle.oracle",
    )
    parser.add_argument(
        "--output-dir",
        default="results/celloracle_hif1a_ko_1500/figures",
    )
    parser.add_argument(
        "--table-dir",
        default="results/celloracle_hif1a_ko_1500/tables",
    )
    parser.add_argument("--target-gene", default="HIF1A")
    parser.add_argument("--max-propagation", type=int, default=5)
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    table_dir = Path(args.table_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    oracle = load_oracle(args.oracle)
    propagation = celloracle_impact_summary(
        oracle,
        args.target_gene,
        max_propagation=args.max_propagation,
        norm_order=1,
    )
    table_path = table_dir / "HIF1A_KO_propagation_mean_delta_x_length_by_cluster.csv"
    propagation.to_csv(table_path, index=False)

    output_base = output_dir / "HIF1A_KO_propagation_magnitude_by_cluster"
    plot_propagation_magnitude(
        propagation,
        "cluster",
        output_base,
        target_cluster="4",
        title="HIF1A KO propagation magnitude by cluster",
        value_column="mean_delta_x_length",
        y_label="Mean delta X length",
    )
    print(output_base.with_suffix(".png").resolve())
    print(output_base.with_suffix(".pdf").resolve())
    print(table_path.resolve())


if __name__ == "__main__":
    main()
