#Usage example: python create_af3_jsoninput.py -i ../../data/benchmark/af3_output/small_db/dataPipelineOut/ -o ../../data/benchmark/benchmark_jsonFiles/small_db
import argparse
from pathlib import Path
import json
import itertools
import copy
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd


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

def standard_json_nomod(data_file, data_nomsa, output_path):
    """ 
    Create 5  combinations of paired, unpaired MSAs and template.
    These consist on taking the standard AF3 output (MSA paired + unpaired) and template selection on
    MSA and performing a permutation of all possible combinations into different json files

    """
    options = ['no', '']

    for paired, unpaired, template in itertools.product(options, repeat=3):
        if paired == "" and unpaired == "" and template == "":
            continue
        if paired == "no" and unpaired == "" and template == "":
            continue
        if paired == "" and unpaired == "no" and template == "":
            continue

        combo = prepare_base(data_file)
        name = combo["name"]

        for s, sD in zip(combo.get("sequences"), data_nomsa.get("sequences")):
            chain = s.get("protein")
            chainD = sD.get("protein")

            if paired == "no":
                chain["pairedMsa"] = ""
            if unpaired == "no":
                chain["unpairedMsa"] = ""
            if template == "no":
                chain["templates"] = []

            if paired == "no" and unpaired == "no" and template=="":
                chain["templates"] = chainD["templates"]

        out_dir = ( Path(output_path) / f"json_{unpaired}unpaired_{paired}paired_{template}template" / name)
        write_af_json(combo, out_dir)

def unpairedWithUniprot(data_uniprotUnpaired = None, output_path = "."):
    """ 
    Create the combination of the unpaired with Uniprot data with notemplate
    """

    template = "no"
    combo = prepare_base(data_uniprotUnpaired)
    name = combo["name"]
    for s in combo.get("sequences"):
        chain = s.get("protein")
        chain["templates"] = []

    out_dir = ( Path(output_path) / f"json_unpairedWithUniprot_nopaired_{template}template" / name)
    write_af_json(combo, out_dir)

def templateOnQuery(data_templateOQ = None, data_standard = None,  output_path = "."):
    """ 
    Create the 4 combinations of the unpaired and paired data given that template was selected on query
    """

    options = ['no', '']

    for paired, unpaired in itertools.product(options, repeat=2):
        combo = prepare_base(data_templateOQ)
        name = combo["name"]
        for sTOQ, s in zip(combo.get("sequences"), data_standard.get("sequences")):
            chain = s.get("protein")
            chainTOQ = sTOQ.get("protein")
            if paired == "":
                chainTOQ["pairedMsa"] = chain["pairedMsa"]
            elif paired == "no":
                chainTOQ["pairedMsa"] = ''
            if unpaired == "":
                chainTOQ["unpairedMsa"] = chain["unpairedMsa"]
            elif unpaired == "no":
                chainTOQ["unpairedMsa"] = ''

        out_dir = (Path(output_path) / f"json_{unpaired}unpaired_{paired}paired_onquerytemplate" / name)
        write_af_json(combo, out_dir)


def reconstruct_protein(mhc = None, tra = None, trb = None, row = None, output_path ='.', folder_name="NAMEFOLDER", description="generated on ..."):

    id_list = ["D", "E", "C", "A"]
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
    jsonf["sequences"] = proteins
    jsonf["modelSeeds"] = [1]
    jsonf["dialect"] = "alphafold3"
    jsonf["version"] = 1

    out_dir = (output_path / folder_name / row["name"])
    write_af_json(jsonf, out_dir)


def process_protein(args):
    protein_name, folderA, folderB, folderC, row, output_path = args
    mhc_name = row["MHCA_aa_id"]
    tra_name = row["TRA_aa_id"]
    trb_name = row["TRB_aa_id"]
    description = f"{mhc_name}{tra_name}{trb_name}"

    mhc_file_A = folderA / mhc_name / f"{mhc_name}_data.json"
    mhc_file_B = folderB / mhc_name / f"{mhc_name}_data.json"
    mhc_file_C = folderC / mhc_name / f"{mhc_name}_data.json"
    tra_file_A = folderA / tra_name / f"{tra_name}_data.json"
    tra_file_B = folderB / tra_name / f"{tra_name}_data.json"
    tra_file_C = folderC / tra_name / f"{tra_name}_data.json"
    trb_file_A = folderA / trb_name / f"{trb_name}_data.json"
    trb_file_B = folderB / trb_name / f"{trb_name}_data.json"
    trb_file_C = folderC / trb_name / f"{trb_name}_data.json"
    
    if not (mhc_file_A.exists() and mhc_file_B.exists() and mhc_file_C.exists() and tra_file_A.exists() and tra_file_B.exists() and tra_file_C.exists() and trb_file_A.exists() and trb_file_B.exists() and trb_file_C.exists()):
        return f"Skipped {protein_name}"
    
    with open(mhc_file_A) as mhcA, open(mhc_file_B) as mhcB, open(mhc_file_C) as mhcC, open(tra_file_A) as traA, open(tra_file_B) as traB, open(tra_file_C) as traC, open(trb_file_A) as trbA, open(trb_file_B) as trbB, open(trb_file_C) as trbC:
        mhcA = json.load(mhcA)
        mhcB = json.load(mhcB)
        mhcC = json.load(mhcC)  
        traA = json.load(traA) 
        traB = json.load(traB)
        traC = json.load(traC)
        trbA = json.load(trbA) 
        trbB = json.load(trbB)  
        trbC = json.load(trbC)  

    reconstruct_protein(mhc = mhcA, tra = traA, trb=trbA, row = row, output_path = output_path, folder_name = "json_unpaired_paired_template", description = description)
    reconstruct_protein(mhc = mhcB, tra = traB, trb=trbB, row = row, output_path = output_path, folder_name = "json_unpairedWithUniprot_nopaired_template", description = description) 
    reconstruct_protein(mhc = mhcC, tra = traC, trb=trbC, row = row, output_path = output_path, folder_name = "json_unpairedWithUniprot_nopaired_onquerytemplate", description = description) 


    reconstructed_file_A = output_path / "json_unpaired_paired_template" / row["name"] / "alphafold_input.json"
    reconstructed_file_B = output_path / "json_unpairedWithUniprot_nopaired_template" / row["name"] / "alphafold_input.json"
    reconstructed_file_C = output_path / "json_unpairedWithUniprot_nopaired_onquerytemplate" / row["name"] / "alphafold_input.json"
    
    with open(reconstructed_file_A) as fA, open(reconstructed_file_B) as fB, open(reconstructed_file_C) as fC:
        fileA = json.load(fA)
        fileB = json.load(fB)
        fileC = json.load(fC)

    standard_json_nomod(fileA, fileC, output_path)
    unpairedWithUniprot(fileC, output_path)
    templateOnQuery(fileC, fileA, output_path)

def main():
    args = parse_args()
    input_dir = Path(args.input)
    og_data = pd.read_csv(args.data)
    folderA = input_dir / "uniprotOn_paired_template_standard"
    folderB = input_dir / "uniprotOn_unpaired_template_standard"
    folderC = input_dir / "uniprotOn_unpaired_template_onquery"

    proteins = og_data["name"]
    tasks = [
        (p, folderA, folderB, folderC, og_data[og_data["name"] == p].iloc[0], args.output)
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
