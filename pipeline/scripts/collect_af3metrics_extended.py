import argparse
from pathlib import Path
import pandas as pd
import json
import numpy as np
import warnings
from concurrent.futures import ProcessPoolExecutor
import os

import biotite.structure.io as strucio

from scoring_utils import get_plddt_score, mean_min_pae_perinterface, mean_plddt_perinterface

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the folder with af3 output containing confidence metrics",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        type=Path,
        required=True,
        help="Suffix for the output_files",
    )
    return parser.parse_args()

def process_complex(args):
    c, suffix = args
    pdb_id = c.name
    ranking_scores = pd.read_csv(c / f"{pdb_id}_ranking_scores.csv")

    # Build fast lookup once
    ranking_lookup = {
        (row.seed, row.sample): float(row.ranking_score)
        for row in ranking_scores.itertuples()
    }

    rows = []
    models = [m for m in c.iterdir() if m.is_dir()]

    for m in models:
        model_name = m.name
        name_split = model_name.split("_")
        seed = int(name_split[0].split("-")[-1])
        sample = int(name_split[1].split("-")[-1])

        af_confidence = ranking_lookup[(seed, sample)]
        dockq_path = m / f"dockQ_metrics_{pdb_id}_{model_name}{suffix}.json"

        if dockq_path.exists():
            with open(dockq_path) as f:
                d = json.load(f)
                dockq = d["DockQ"]
        else:
            print(f"File {dockq_path} does not exist.")
            dockq = None 

        with open(m / f"{pdb_id}_{model_name}_summary_confidences.json") as f:
            afsumm = json.load(f)

        with open(m / f"{pdb_id}_{model_name}_confidences.json") as f:
            af = json.load(f)
            
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            structure = strucio.load_structure(
                m / f"{pdb_id}_{model_name}_model.cif"
            )
            structure.set_annotation("atom_id", np.arange(1, len(structure)+1))

        df_plddt = mean_plddt_perinterface(structure, af)
        df_min_pae, df_mean_pae = mean_min_pae_perinterface(structure, af)
        plddt_cdrpep = get_plddt_score(structure, plddt_scores=af["atom_plddts"])

        rows.append({
            "pdb_id": pdb_id,
            "seed": seed,
            "sample": sample,
            "sample_name": model_name,
            "af_confidence": af_confidence,
            "chain_iptm": afsumm["chain_iptm"],
            "chain_pair_iptm": afsumm["chain_pair_iptm"],
            "chain_pair_pae_min": afsumm["chain_pair_pae_min"],
            "chain_ptm": afsumm["chain_ptm"],
            "iptm": float(afsumm["iptm"]),
            "ptm": float(afsumm["ptm"]),
            "mean_plddt": float(pd.Series(af["atom_plddts"]).mean()),
            "mean_plddt_perinterface": df_plddt.to_dict(),
            "min_pae_perinterface": df_min_pae.to_dict(),
            "mean_pae_perinterface": df_mean_pae.to_dict(),
            "plddt_cdrpep": plddt_cdrpep,
            "dockq":dockq
        })
    df_pdb = pd.DataFrame(rows)
    df_pdb.to_csv(c / f"af3_extended_metrics{suffix}.csv", index=False)

    return df_pdb

def main():
    args = parse_args()
    input_path = args.input
    suffix = args.suffix

    complexes = [c for c in input_path.iterdir() if c.is_dir()]
    tasks = [(c, suffix) for c in complexes]
    dfs = []
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        for df in executor.map(process_complex, tasks):
            dfs.append(df)

    df_metrics = pd.concat(dfs, ignore_index=True)
    df_metrics.to_csv(input_path / f"collected_af3metrics{suffix}.csv", index=False)

if __name__ == "__main__":
    main()