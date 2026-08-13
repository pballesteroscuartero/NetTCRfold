#!/usr/bin/env python3
#Usage mode: python collect_rerankingMetrics.py -i ../../data/nettcrstruc/rerank_input/inputModels/ -o ../../data/af3_output/small_db/structInference/

import pandas as pd
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the folder with gnn rescoring",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output folder where to save gnn rescoring",
    )
    parser.add_argument(
        "-s",
        '--suffix',
        type=str,
        default="",
        help="Suffix to add to the output file name (default: '')"
    )
    return parser.parse_args()


def load_and_merge_gnn_reranking(args):
    complex_dir, param_comb = args
    pdb_id = complex_dir.name
    df_if1 = complex_dir / "rescore_ensemble_gvp_if1.csv"
    df_gvp = complex_dir / "rescore_ensemble_gvp.csv"

    if not df_if1.exists() or not df_gvp.exists():
        print(f"Missing files for {pdb_id} in {complex_dir}. Skipping.")
        return None

    df_if1 = pd.read_csv(df_if1)
    df_gvp = pd.read_csv(df_gvp)

    for df in (df_if1, df_gvp):
        df["param_comb"] = param_comb
        df["pdb_id"] = pdb_id
        df["sample_name"] = (
            df["name"].str.split("_").str[-2]
            + "_"
            + df["name"].str.split("_").str[-1]
        )

    df_if1 = df_if1[["pdb_id", "sample_name", "param_comb", "gnn_ens"]]
    df_gvp = df_gvp[["pdb_id", "sample_name", "param_comb", "gnn_ens"]]

    df_nettcrstruct = pd.merge(
        df_if1,
        df_gvp,
        on=["pdb_id", "sample_name", "param_comb"],
        suffixes=("_if1", "")
    )

    return df_nettcrstruct

def main():
    args = parse_args()
    input_path = args.input
    output_path = args.output
    
    folders = [f for f in input_path.iterdir() if f.is_dir() and "json_" in(f.name)]
    tasks = []

    for folder in folders:
        job_splits = [js for js in folder.iterdir() if js.is_dir() and "job_split" in js.name]
        for js in job_splits:
            complexes = [c for c in js.iterdir() if c.is_dir()]
            tasks.extend([(c, folder.name) for c in complexes])
        #complexes = [c for c in folder.iterdir() if c.is_dir()]
        #tasks.extend([(c, folder.name) for c in complexes])
    with ProcessPoolExecutor(max_workers=10) as executor:
        dfs = list(executor.map(load_and_merge_gnn_reranking, tasks))

    dfs = [df for df in dfs if df is not None]
    df_rerank = pd.concat(dfs, ignore_index=True)
    out_file = output_path / f"all_gnnreranked{args.suffix}.csv"
    df_rerank.to_csv(out_file, index=False)

if __name__ == "__main__":
    main()

