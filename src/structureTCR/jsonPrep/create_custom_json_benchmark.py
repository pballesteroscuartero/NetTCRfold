#Generates the full MSA/template combination sweep used for benchmarking (paired/unpaired MSA x template mode).
#Usage example: python create_custom_json_benchmark.py -i ../../data/benchmark/af3_output/small_db/dataPipelineOut/ -o ../../data/benchmark/benchmark_jsonFiles/small_db

import argparse
import pandas as pd
from pathlib import Path
import json
import itertools
import copy
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


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

def templateOnQuery(file, output_path = "."):
    
    options = ['no', '']
    for paired, unpaired in itertools.product(options, repeat=2):
        if paired == "" and unpaired == "":
            continue
        if unpaired == "no" and paired == "":
            continue
        combo = prepare_base(file)
        name = combo["name"]
        for s in combo.get("sequences"):
            chain = s.get("protein")
            if paired == "no":
                chain["pairedMsa"] = ""
            if unpaired == "no":
                chain["unpairedMsa"] = ""

        out_dir = (Path(output_path) / f"json_{unpaired}unpaired_{paired}paired_onquerytemplate" / name)
        write_af_json(combo, out_dir)
        
def reconstruct_protein(mhc = None, tra = None, trb = None, mhcb = None, row = None, output_path ='.', folder_name="NAMEFOLDER", description="generated on ..."):

    id_list = ["D", "E", "C", "A"]
    if mhcb is not None:
        id_list.append("B")
    jsonf = {}
    jsonf["name"] = row["name"]
    proteins = []
    for id in id_list:
        if id == "D": #This is the TRA
            sequences = tra.get("sequences")
            if len(sequences) != 1:
                raise ValueError(
                    f"{id}: expected exactly 1 sequence, found {len(sequences)}"
                )
            for s in sequences:
                protein_info = s.get("protein")
            protein = {
                    "protein": protein_info
                }
            proteins.append(protein)
        
        elif id == "E": #This is the TRB
            sequences = trb.get("sequences")
            if len(sequences) != 1:
                raise ValueError(
                    f"{id}: expected exactly 1 sequence, found {len(sequences)}"
                )
            for s in sequences:
                protein_info = s.get("protein")
            protein = {
                    "protein": protein_info
                }
            proteins.append(protein)

        elif id == "C": #This is peptide
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
        elif id == "A": #This is mhc
            sequences = mhc.get("sequences")
            if len(sequences) != 1:
                raise ValueError(
                    f"{id}: expected exactly 1 sequence, found {len(sequences)}"
                )
            for s in sequences:
                protein_info = s.get("protein")
            protein = {
                    "protein": protein_info
                }
            proteins.append(protein)
        elif id == "B": #This is mhcb
            sequences = mhcb.get("sequences")
            if len(sequences) != 1:
                raise ValueError(
                    f"{id}: expected exactly 1 sequence, found {len(sequences)}"
                )
            for s in sequences:
                protein_info = s.get("protein")
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
    protein_name, folderA, folderB, row, output_path = args
    mhc_name = row["MHCA_aa_id"]
    tra_name = row["TRA_aa_id"]
    trb_name = row["TRB_aa_id"]

    mhcb_name = row.get("MHCB_aa_id")
    has_mhcb = mhcb_name is not None and pd.notna(mhcb_name)
    description = f"{mhc_name}{tra_name}{trb_name}"
    
    mhc_file_A = folderA / mhc_name / f"{mhc_name}_data.json"
    mhc_file_B = folderB / mhc_name / f"{mhc_name}_data.json"
    tra_file_A = folderA / tra_name / f"{tra_name}_data.json"
    tra_file_B = folderB / tra_name / f"{tra_name}_data.json"
    trb_file_A = folderA / trb_name / f"{trb_name}_data.json"
    trb_file_B = folderB / trb_name / f"{trb_name}_data.json"
    
    if has_mhcb:
        description += f"{mhcb_name}"
        mhcb_file_A = folderA / mhcb_name / f"{mhcb_name}_data.json"
        mhcb_file_B = folderB / mhcb_name / f"{mhcb_name}_data.json"

    if not (mhc_file_A.exists() and mhc_file_B.exists() and tra_file_A.exists() and tra_file_B.exists() and trb_file_A.exists() and trb_file_B.exists()):
        if has_mhcb:
            if not (mhcb_file_A.exists() and mhcb_file_B.exists()):
                return f"Skipped {protein_name} because of missing MHCb file"
        else:
            return f"Skipped {protein_name}"
    mhcbA = None
    mhcbB = None
    with open(mhc_file_A) as mhcA, open(mhc_file_B) as mhcB, open(tra_file_A) as traA, open(tra_file_B) as traB, open(trb_file_A) as trbA, open(trb_file_B) as trbB:
        if has_mhcb:
            with open(mhcb_file_A) as mhcbA, open(mhcb_file_B) as mhcbB:
                mhcbA = json.load(mhcbA)
                mhcbB = json.load(mhcbB)
        mhcA = json.load(mhcA)
        mhcB = json.load(mhcB) 
        traA = json.load(traA) 
        traB = json.load(traB)
        trbA = json.load(trbA) 
        trbB = json.load(trbB)  

    reconstruct_protein(mhc = mhcA, tra = traA, trb=trbA, mhcb= mhcbA, row = row, output_path = output_path, folder_name = "json_unpaired_paired_onquerytemplate", description = description)
    reconstruct_protein(mhc = mhcB, tra = traB, trb=trbB, mhcb= mhcbB, row = row, output_path = output_path, folder_name = "json_unpairedWithUniprot_nopaired_onquerytemplate", description = description)

    reconstructed_file = output_path / "json_unpaired_paired_onquerytemplate" / row["name"] / "alphafold_input.json"
    with open(reconstructed_file) as f:
        file = json.load(f)
    templateOnQuery(file, output_path) 


def main():
    args = parse_args()
    input_dir = Path(args.input)
    og_data = pd.read_csv(args.data)
    folderA = input_dir / "uniprotOn_paired_template_onquery"
    folderB = input_dir / "uniprotOn_unpaired_template_onquery"

    proteins = og_data["name"]
    tasks = [
        (p, folderA, folderB, og_data[og_data["name"] == p].iloc[0], args.output)
        for p in proteins
    ]
    with ProcessPoolExecutor(max_workers=20) as ex:
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
