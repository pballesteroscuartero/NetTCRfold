#Merges the per-split collected_af3metrics CSVs produced by collect_af3metrics_extended_parallel.py into one file.
import argparse
import pandas as pd
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Folder containing the collected_af3metrics_split*.csv files",
    )
    parser.add_argument(
        "-n",
        "--num-splits",
        type=int,
        required=True,
        help="Number of splits to merge (matches --num_splits used when collecting metrics)",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        default="",
        help="Suffix used when the split files were generated (e.g. _dockQonpMHC)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = args.input
    suffix = args.suffix

    dfs = [pd.read_csv(input_path / f"collected_af3metrics_split{i}{suffix}.csv") for i in range(args.num_splits)]
    df_all = pd.concat(dfs, ignore_index=True)
    df_all.to_csv(input_path / f"collected_af3metrics{suffix}.csv", index=False)


if __name__ == "__main__":
    main()