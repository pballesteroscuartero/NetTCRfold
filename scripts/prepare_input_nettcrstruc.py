import shutil
import pandas as pd
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial

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
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the folder where nettcrstruct formatted input will be stored",
    )
    return parser.parse_args()


def process_complex(complex_dir, structcrinput_dir):
    if not complex_dir.is_dir():
        return f"Skipped {complex_dir.name}"

    complex_name = complex_dir.name
    #print(f"Processing complex {complex_name}")

    run_input_dir = structcrinput_dir / complex_name
    run_input_dir.mkdir(parents=True, exist_ok=True)

    seed_dirs = [d for d in complex_dir.iterdir() if d.is_dir() and "seed" in d.name]

    pdb_records = []
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

    ranking_csv = complex_dir / f"{complex_name}_ranking_scores.csv"
    if not ranking_csv.exists():
        print(f"WARNING: Ranking score file missing for {complex_name}")
        return f"Skipped {complex_name}"

    df = pd.read_csv(ranking_csv)
    df["name"] = df.apply(
        lambda x: f"{complex_name}_seed-{int(x['seed'])}_sample-{int(x['sample'])}", axis=1
    )
    df.rename(columns={"ranking_score": "confidence"}, inplace=True)
    df2 = df[["name", "confidence"]]
    df2.to_csv(run_input_dir / "model_scores.txt", sep="\t", index=False)

    #print(f"Prepared {complex_name}: {len(pdb_records)} models")
    return complex_name


def main():
    args = parse_args()
    input_folder = Path(args.input)
    output_folder = Path(args.output)
    folder_name = input_folder.name

    structcrinput_dir = output_folder / folder_name
    structcrinput_dir.mkdir(parents=True, exist_ok=True)

    complex_dirs = [f for f in input_folder.iterdir() if f.is_dir()]
    
    process_func = partial(process_complex, structcrinput_dir=structcrinput_dir)
    with ProcessPoolExecutor(max_workers=20) as executor:
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