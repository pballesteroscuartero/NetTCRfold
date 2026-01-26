#Usage example: python create_af3_jsoninput.py -i ../../data/benchmark/af3_output/small_db/dataPipelineOut/ -o ../../data/benchmark/benchmark_jsonFiles/small_db
import argparse
from pathlib import Path
import json
import itertools
import copy
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor



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

def standard_json_nomod(pdb_id, data_file, output_path):
    """ 
    Create the 8 combinations of paired, unpaired MSAs and template.
    These consist on taking the standard AF3 output (MSA paired + unpaired) and template selection on
    MSA and performing a permutation of all possible combinations into different json files

    """
    options = ['no', '']

    for paired, unpaired, template in itertools.product(options, repeat=3):
        combo = prepare_base(data_file)
        for s in combo.get("sequences"):
            chain = s.get("protein")
            if paired == "no":
                chain["pairedMsa"] = ""
            if unpaired == "no":
                chain["unpairedMsa"] = ""
            if template == "no":
                chain["templates"] = []

        out_dir = (
            output_path /
            f"json_{unpaired}unpaired_{paired}paired_{template}template" /
            pdb_id
        )
        write_af_json(combo, out_dir)

def unpairedWithUniprot(pdb_id, data_uniprotUnpaired = None, data_templateOQ = None, output_path = "."):
    """ 
    Create the 3 combinations of the unpaired with Uniprot data and the 3 template possibilities (yes, no, templateOnQuery)
    """

    template_options = ['no', '', 'onquery']

    for template in template_options:
        combo = prepare_base(data_uniprotUnpaired)
        for s, sTOQ in zip(combo.get("sequences"), data_templateOQ.get("sequences")):
            chain = s.get("protein")
            chainTOQ = sTOQ.get("protein")
            if template == "no":
                chain["templates"] = []
            elif template == "onquery":
                chain["templates"] = chainTOQ["templates"]

        out_dir = (
            output_path /
            f"json_unpairedWithUniprot_nopaired_{template}template" /
            pdb_id
        )
        write_af_json(combo, out_dir)

def templateOnQuery(pdb_id, data_templateOQ = None, data_standard = None,  output_path = "."):
    """ 
    Create the 4 combinations of the unpaired and paired data given that template was selected on query
    """

    options = ['no', '']

    for paired, unpaired in itertools.product(options, repeat=2):
        combo = prepare_base(data_templateOQ)

        for sTOQ, s in zip(combo.get("sequences"), data_standard.get("sequences")):
            chain = s.get("protein")
            chainTOQ = sTOQ.get("protein")
            if paired == "":
                chainTOQ["pairedMsa"] = chain["pairedMsa"]
            if unpaired == "":
                chainTOQ["unpairedMsa"] = chain["unpairedMsa"]

        out_dir = (
            output_path /
            f"json_{unpaired}unpaired_{paired}paired_onquerytemplate" /
            pdb_id
        )
        write_af_json(combo, out_dir)

def process_protein(args):
    protein_name, folderA, folderB, folderC, output_path = args

    fileA = folderA / protein_name / f"{protein_name}_data.json"
    fileB = folderB / protein_name / f"{protein_name}_data.json"
    fileC = folderC / protein_name / f"{protein_name}_data.json"

    if not (fileA.exists() and fileB.exists() and fileC.exists()):
        return f"Skipped {protein_name}"

    with open(fileA) as fA, open(fileB) as fB, open(fileC) as fC:
        dataA = json.load(fA)
        dataB = json.load(fB)
        dataC = json.load(fC)

    standard_json_nomod(protein_name, dataA, output_path)
    unpairedWithUniprot(protein_name, dataB, dataC, output_path)
    templateOnQuery(protein_name, dataC, dataA, output_path)

    return protein_name


def main():
    args = parse_args()
    input_dir = Path(args.input)

    folderA = input_dir / "uniprotOn_paired_template_standard"
    folderB = input_dir / "uniprotOn_unpaired_template_standard"
    folderC = input_dir / "uniprotOn_unpaired_template_onquery"

    proteins = [p.name for p in folderA.iterdir() if p.is_dir()]

    tasks = [
        (p, folderA, folderB, folderC, args.output)
        for p in proteins
    ]

    with ProcessPoolExecutor() as ex:
        ex.map(process_protein, tasks)
    

if __name__ == "__main__":
    main()
