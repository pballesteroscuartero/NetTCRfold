import shutil
import pandas as pd
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import ast
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the folder containing the param_comb folder with the different complexes ",
    )
    parser.add_argument(
        "-n",
        "--name",
        required=True,
        help="Folder name",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the folder where nettcrstruct formatted input will be stored",
    )
    return parser.parse_args()


def flatten_dict(d, key_map, mode="offdiag"):

    if isinstance(d, str):
        try:
            d = ast.literal_eval(d)
        except Exception:
            d = eval(d.replace('nan', 'np.nan'))
    keys_in_row = [k for k in key_map.keys() if k in d]
    values = []
    key_pairs = []
    for i, r_key in enumerate(keys_in_row):
        for j, c_key in enumerate(keys_in_row):
            if mode =="upper":
                if j > i:
                    values.append(d[r_key][c_key])
                    key_pairs.append((key_map[r_key], key_map[c_key]))
            elif mode =="offdiag":
                if i != j:
                    values.append(d[r_key][c_key])
                    key_pairs.append((key_map[r_key], key_map[c_key]))

    return values, key_pairs

def extract_dictionary_into_df(df, row, idx, metrics, dict_key_map, mode_flatten = "offdiag"):
    for metric in metrics:
        matrix = row[metric]
        vals, order =flatten_dict(matrix, dict_key_map, mode=mode_flatten)
        for col_name, val in zip(order, vals):
            colname = f"{col_name[0]}_{col_name[1]}_{metric}"
            df.at[idx, colname] = val

def process_complex(complex_dir, structcrinput_dir):
    if not complex_dir.is_dir():
        return f"Skipped {complex_dir.name}"

    complex_name = complex_dir.name
    #print(f"Processing complex {complex_name}")

    run_input_dir = structcrinput_dir / complex_name
    run_input_dir.mkdir(parents=True, exist_ok=True)

    seed_dirs = [d for d in complex_dir.iterdir() if d.is_dir() and "seed" in d.name]

    pdb_records = []
    ipsae_paths = []

    for sd in seed_dirs:
        pdb_files = list(sd.glob("*model.cif"))
        if len(pdb_files) != 1:
            print(f"WARNING: {sd} does not have exactly one PDB file.")
            continue

        pdb_file = pdb_files[0]
        model_name = pdb_file.stem.replace("_model", "")
        new_pdb_name = f"{model_name}.cif"

        shutil.copy(pdb_file, run_input_dir / new_pdb_name)
        pdb_records.append(model_name)

    # ranking_csv = complex_dir / f"af3_extended_metrics_dockQonpMHC.csv"
    # if not ranking_csv.exists():
    #     print(f"WARNING: Ranking score file missing for {complex_name}")
    #     return f"Skipped {complex_name}"
    
    # df = pd.read_csv(ranking_csv)

    # df_subset = df[["pdb_id", "sample_name", "af_confidence", "cdr_metric_mean_chain", "ipsae_d0dom"]]
    # metrics = ["cdr_metric_mean_chain", "ipsae_d0dom"]
    # dict_key_map = {'A': 'MHCA', 'C': 'pep', 'D': 'TRA', 'E': 'TRB'}
    # for idx, row in df_subset.iterrows():
    #     extract_dictionary_into_df(df_subset, row, idx, metrics, dict_key_map, mode_flatten="offdiag")
    # df_subset["TCR_pep_ipsae"] = df_subset[["TRB_pep_ipsae_d0dom", "TRA_pep_ipsae_d0dom", "pep_TRB_ipsae_d0dom", "pep_TRA_ipsae_d0dom"]].max(axis=1)
    # df_subset["MHCA_CDR12_ipsae"] = df_subset[["MHCA_TRA_cdr_metric_mean_chain", "MHCA_TRB_cdr_metric_mean_chain", "TRA_MHCA_cdr_metric_mean_chain", "TRB_MHCA_cdr_metric_mean_chain"]].max(axis=1)
    # df_subset["name"] = df_subset["pdb_id"] + "_" + df_subset["sample_name"]
    # cols_keep = ["name","af_confidence", "TCR_pep_ipsae", "MHCA_CDR12_ipsae"]
    # df_final = df_subset[cols_keep]
    # SCALE = 10  # PAE of 10Å gives score of 0.5 UNDERSTAND BETTER
    # df["MHCA_CDR12_ipsae"] = 1 / (1 + df["MHCA_CDR12_ipsae"] / SCALE)
    # epsilon = 0.01
    # df["TCR_pep_ipsae"] = [v if v > 0 else epsilon for v in df["TCR_pep_ipsae"]]
    # df["MHCA_CDR12_ipsae"] = [v if v > 0 else epsilon for v in df["MHCA_CDR12_ipsae"]]

    # df_final.rename(columns={"af_confidence": "confidence"}, inplace=True)
    # #df2 = df[["name", "confidence"]]
    # df_final.to_csv(run_input_dir / "model_scores.txt", sep="\t", index=False)
    return complex_name


def main():
    args = parse_args()
    input_folder = Path(args.input)
    output_folder = Path(args.output)
    folder_name = args.name

    structcrinput_dir = output_folder / folder_name
    structcrinput_dir.mkdir(parents=True, exist_ok=True)

    complex_dirs = [f for f in input_folder.iterdir() if f.is_dir()]
    
    process_func = partial(process_complex, structcrinput_dir=structcrinput_dir)
    with ProcessPoolExecutor(max_workers=5) as executor:
        executor.map(process_func, complex_dirs)
    print("All complexes processed")


# def parse_args() -> argparse.Namespace:
#     parser = argparse.ArgumentParser()
#     parser.add_argument(
#         "-i",
#         "--input",
#         type=Path,
#         required=True,
#         help="Path to the folder containing the param_comb folder with the different complexes ",
#     )
#     return parser.parse_args()

# def main():
#     args = parse_args()
#     input_folder = args.input
#     folder_name = input_folder.name
#     structcrinput_dir = Path(f"/home/projects2/pbacu/projects/structureTCR/schumacherDataset/data/nettcrstruc/rerank_input/inputModels/{folder_name}")
#     structcrinput_dir.mkdir(parents=True, exist_ok=True)
    
#     for complex_dir in input_folder.iterdir():
#         if not complex_dir.is_dir():
#             continue

#         complex_name = complex_dir.name
#         print(f"Processing complex {complex_name}")

#         run_input_dir = structcrinput_dir / complex_name
#         run_input_dir.mkdir(exist_ok=True)
        
#         seed_dirs = [d for d in complex_dir.iterdir() if d.is_dir() and "seed" in d.name]

#         pdb_records = []
#         for sd in seed_dirs:
#             pdb_files = list(sd.glob("*model.cif"))
#             if len(pdb_files) != 1:
#                 print(f"WARNING: {sd} does not have exactly one PDB file.")
#                 continue

#             pdb_file = pdb_files[0]
#             model_name = pdb_file.stem.replace("_model", "")  # remove "_model"
#             new_pdb_name = f"{model_name}.cif"

#             shutil.copy(pdb_file, run_input_dir / new_pdb_name)
#             pdb_records.append(model_name)

#         ranking_csv = complex_dir / f"{complex_name}_ranking_scores.csv"
#         if not ranking_csv.exists():
#             print(f"WARNING: Ranking score file missing for {complex_name}")
#             continue

#         df = pd.read_csv(ranking_csv)
#         df["name"] = df.apply(lambda x: f"{complex_name}_seed-{int(x['seed'])}_sample-{int(x['sample'])}", axis=1)
#         df.rename(columns={"ranking_score": "confidence"}, inplace=True)

#         # Keep only relevant fields
#         df2 = df[["name", "confidence"]]

#         # Save as model_scores.txt
#         df2.to_csv(run_input_dir / "model_scores.txt", sep="\t", index=False)

#         print(f"Prepared {complex_name}: {len(pdb_records)} models")

#     print("All complexes processed.")

if __name__ == "__main__":
    main()