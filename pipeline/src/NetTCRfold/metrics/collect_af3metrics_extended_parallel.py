import argparse
from pathlib import Path
import pandas as pd
import json
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"


import numpy as np
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess

from tqdm import tqdm
import biotite.structure.io as strucio

from NetTCRfold.metrics.scoring_utils import is_interface, pae_metrics_perchain_pair, parse_ipsae, get_plddt_score
from NetTCRfold.metrics.structure_utils import extract_cdrs_from_structure
from NetTCRfold.metrics.ipsae_function import compute_ipsae

# Must match the -s suffix workflow.sh passes to computeDockq.py.
DOCKQ_SUFFIX = "_dockQonpMHC"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-s", "--suffix", required=True)
    parser.add_argument(
        "--split_idx",
        type=int,
        default=0,
        help="Index of this split (0-based)",
    )
    parser.add_argument(
        "--num_splits",
        type=int,
        default=1,
        help="Total number of splits",
    )
    return parser.parse_args()

def has_any_dockq(complexes):
    """Whether any DockQ output exists yet for this folder (short-circuits on the first hit)."""
    for c in complexes:
        pdb_id = c.name
        for m in c.iterdir():
            if m.is_dir() and (m / f"dockQ_metrics_{pdb_id}_{m.name}{DOCKQ_SUFFIX}.json").exists():
                return True
    return False


def process_complex(args):
    c, has_dockq = args
    pdb_id = c.name
    ranking_scores = pd.read_csv(c / f"{pdb_id}_ranking_scores.csv")

    # Build fast lookup once
    ranking_lookup = {
        (row.seed, row.sample): float(row.ranking_score)
        for row in ranking_scores.itertuples()
    }

    rows = []
    models = [m for m in c.iterdir() if m.is_dir()]
    cdrs_cache = None

    for m in models:
        model_name = m.name
        name_split = model_name.split("_")
        seed = int(name_split[0].split("-")[-1])
        sample = int(name_split[1].split("-")[-1])

        af_confidence = ranking_lookup[(seed, sample)]
        model_summary_confidences_path = m / f"{pdb_id}_{model_name}_summary_confidences.json"
        model_confidences_path = m / f"{pdb_id}_{model_name}_confidences.json"
        structure_path = m / f"{pdb_id}_{model_name}_model.cif"

        dockq = None
        if has_dockq:
            dockq_path = m / f"dockQ_metrics_{pdb_id}_{model_name}{DOCKQ_SUFFIX}.json"
            if dockq_path.exists():
                with open(dockq_path) as f:
                    dockq = json.load(f)["DockQ"]
            else:
                print(f"No dockQ file found for {pdb_id}/{model_name}.")

        with open(model_summary_confidences_path) as f:
            afsumm = json.load(f)

        with open(model_confidences_path) as f:
            af = json.load(f)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            structure = strucio.load_structure( structure_path)
            structure.set_annotation("atom_id", np.arange(1, len(structure)+1))
            structure_noh = structure[structure.element != 'H']

            noh_mask = structure.element != 'H'
            af_arrays = {
                "atom_plddts": np.array(af["atom_plddts"])[noh_mask],
                "token_chain_ids": np.array(af["token_chain_ids"]),
                "token_res_ids": np.array(af["token_res_ids"]),
                "pae": np.array(af["pae"]),
            }

        df_ipsae = compute_ipsae(af, afsumm, structure, structure_path, 10, 10)
   
        if cdrs_cache is None:
            cdrs_cache = extract_cdrs_from_structure(structure_noh)
            
        df_interface = is_interface(structure_noh)
        plddt_cdrpep = get_plddt_score(structure_noh, plddt_scores=af_arrays["atom_plddts"]) #Corrected with the way of extracting CDRs
        df_cdr_metric_mean = pae_metrics_perchain_pair(structure_noh, af_arrays, cdrs = cdrs_cache)
        metric_dict_ipsae = parse_ipsae(df_ipsae)
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
            "is_interface": df_interface.to_dict(),
            "plddt_cdrpep": plddt_cdrpep,
            "dockq":dockq,
            "cdr_metric_mean_chain": df_cdr_metric_mean.to_dict(),
            "ipsae": metric_dict_ipsae["ipSAE"],
            "ipsae_d0chn": metric_dict_ipsae["ipSAE_d0chn"],
            "ipsae_d0dom": metric_dict_ipsae["ipSAE_d0dom"],            
        })

    df_pdb = pd.DataFrame(rows)
    return df_pdb


def main():
    args = parse_args()
    input_path = args.input
    suffix = args.suffix

    complexes = sorted([c for c in input_path.iterdir() if c.is_dir()])  # sort for determinism across jobs

    # Partition into num_splits contiguous chunks, this job handles only split_idx
    complexes = complexes[args.split_idx::args.num_splits]

    has_dockq = has_any_dockq(complexes)

    n_workers = 4
    joint_path = input_path / f"collected_af3metrics_split{args.split_idx}{suffix}.csv"
    first_write = True

    tasks = [(c, has_dockq) for c in complexes]

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(process_complex, task) for task in tasks]
        for future in tqdm(as_completed(futures), total=len(tasks), desc=f"Split {args.split_idx}/{args.num_splits}"):
            df = future.result()
            if df.empty:
                continue
            pdb_id = df["pdb_id"].iloc[0]
            c = input_path / pdb_id
            out_path = c / f"af3_extended_metrics{suffix}.csv"
            df.to_csv(out_path, index=False)
            df.to_csv(joint_path, mode="a", header=first_write, index=False)
            first_write = False


if __name__ == "__main__":
    main()
