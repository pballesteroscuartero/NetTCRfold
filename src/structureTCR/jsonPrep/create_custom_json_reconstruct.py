#Reconstructs the JSON input for the method used in the final publication (onquery template selection, not the full benchmark sweep).
#Usage example: python create_custom_json_reconstruct.py -i ../../data/benchmark/af3_output/small_db/dataPipelineOut/ -o ../../data/benchmark/benchmark_jsonFiles/small_db -d data.csv
import argparse
from pathlib import Path
import json
import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

FOLDER_NAME = "json_unpaired_nopaired_onquerytemplate"
INPUT_SUBFOLDER = "uniprotOn_paired_template_onquery"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the folders of the data pipeline output",
    )
    parser.add_argument(
        "-d",
        "--data",
        type=Path,
        required=True,
        help="Path to csv with original datapoints and unique sequence identifiers",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output file with the results. Provide here the general output folder. Extra folders will be generated for param_combinations",
    )
    return parser.parse_args()


def write_af_json(base_data, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "alphafold_input.json", "w") as wf:
        json.dump(base_data, wf, indent=4)


def prepare_base(data):
    base = copy.deepcopy(data)
    base["modelSeeds"] = [1]
    base["version"] = 1
    return base


def reconstruct_protein(mhc=None, tra=None, trb=None, row=None, output_path='.', folder_name="NAMEFOLDER", description="generated on ..."):

    id_list = ["D", "E", "C", "A"]
    jsonf = {}
    jsonf["name"] = row["name"]
    proteins = []
    for id in id_list:
        if id == "D":  # This is the TRA
            sequences = tra.get("sequences")
            if len(sequences) != 1:
                raise ValueError(
                    f"{id}: expected exactly 1 sequence, found {len(sequences)}"
                )
            for s in sequences:
                protein_info = s.get("protein")
                protein_info["pairedMsa"] = ""
            protein = {
                    "protein": protein_info
                }
            proteins.append(protein)

        elif id == "E":  # This is the TRB
            sequences = trb.get("sequences")
            if len(sequences) != 1:
                raise ValueError(
                    f"{id}: expected exactly 1 sequence, found {len(sequences)}"
                )
            for s in sequences:
                protein_info = s.get("protein")
                protein_info["pairedMsa"] = ""
            protein = {
                    "protein": protein_info
                }
            proteins.append(protein)

        elif id == "C":  # This is peptide
            protein = {
                    "protein": {
                        "id": id,
                        "sequence": row["peptide"],
                        "unpairedMsa": "",
                        "pairedMsa": "",
                        "templates": [],
                        "description": description

                    }
                }
            proteins.append(protein)
        elif id == "A":  # This is mhc
            sequences = mhc.get("sequences")
            if len(sequences) != 1:
                raise ValueError(
                    f"{id}: expected exactly 1 sequence, found {len(sequences)}"
                )
            for s in sequences:
                protein_info = s.get("protein")
                protein_info["pairedMsa"] = ""
            protein = {
                    "protein": protein_info
                }
            proteins.append(protein)
    jsonf["sequences"] = proteins
    jsonf["modelSeeds"] = [1]
    jsonf["dialect"] = "alphafold3"
    jsonf["version"] = 1

    out_dir = (output_path / folder_name / row["name"])
    write_af_json(jsonf, out_dir)


def process_protein(args):
    protein_name, folderA, row, output_path = args
    mhc_name = row["MHCA_aa_id"]
    tra_name = row["TRA_aa_id"]
    trb_name = row["TRB_aa_id"]
    description = f"{mhc_name}{tra_name}{trb_name}"

    mhc_file_A = folderA / mhc_name / f"{mhc_name}_data.json"
    tra_file_A = folderA / tra_name / f"{tra_name}_data.json"
    trb_file_A = folderA / trb_name / f"{trb_name}_data.json"

    if not (mhc_file_A.exists() and tra_file_A.exists() and trb_file_A.exists()):
        return f"Skipped {protein_name}"

    with open(mhc_file_A) as mhcA, open(tra_file_A) as traA, open(trb_file_A) as trbA:
        mhcA = json.load(mhcA)
        traA = json.load(traA)
        trbA = json.load(trbA)

    reconstruct_protein(mhc=mhcA, tra=traA, trb=trbA, row=row, output_path=output_path, folder_name=FOLDER_NAME, description=description)


def main():
    args = parse_args()
    input_dir = Path(args.input)
    og_data = pd.read_csv(args.data)
    folderA = input_dir / INPUT_SUBFOLDER

    proteins = og_data["name"]
    tasks = [
        (p, folderA, og_data[og_data["name"] == p].iloc[0], args.output)
        for p in proteins
    ]
    with ProcessPoolExecutor(max_workers=10) as ex:
        futures = [ex.submit(process_protein, task) for task in tasks]
        for future in as_completed(futures):
            try:
                result = future.result()
                if result is not None:
                    print(result)
            except Exception as e:
                print("Worker failed with exception:")
                raise  # re-raise to crash loudly


if __name__ == "__main__":
    main()
