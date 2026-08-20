#Builds AF3 input JSON from data-generation output for a user-specified set of MSA/template
#combinations. Each combination is "<msaMode>_<templateMode>":
#  msaMode:      unpairedMSA, pairedMSA, fullMSA, noMSA
#  templateMode: onquery, standard, no

#Usage example: python create_custom_json.py -i .../dataPipelineOut -o .../customJSON -d data_hla_withid.csv -c "unpairedMSA_onquery pairedMSA_onquery"
import argparse
from pathlib import Path
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd


MSA_MODES = {"unpaired", "paired", "full", "no"}
TEMPLATE_MODES = {"onquery", "standard", "no"}

#Which MSA fields to null out on the TRA/TRB/MHC chains for a given MSA mode.
MSA_NULL_FIELDS = {
    "unpaired": {"pairedMsa"},
    "paired": {"unpairedMsa"},
    "full": set(),
    "no": {"pairedMsa", "unpairedMsa"},
}

TEMPLATE_INPUT_SUBFOLDER = {
    "onquery": "template_onquery",
    "standard": "template_standard",
}

#MSA types the data-generation step (step 2) can be run with. "full" computes both
#unpairedMsa and pairedMsa; "unpaired"/"paired" skip searching for the other type entirely.
DATAGEN_MSA_TYPES = ("unpaired", "paired", "full")

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the data-generation output (parent of msa_{unpaired,paired,full}_template_{onquery,standard})",
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
        help="Path to the output folder. A subfolder is created per combination (json_<combination>).",
    )
    parser.add_argument(
        "-c",
        "--combinations",
        default="unpaired_onquery",
        help=(
            "Space-separated list of <msaMode>_<templateMode> combinations to generate. "
            f"msaMode: {sorted(MSA_MODES)}. templateMode: {sorted(TEMPLATE_MODES)}. "
            "Example: 'unpaired_onquery paired_pairedMSA_onquery full_no'"
        ),
    )
    return parser.parse_args()


def parse_combination(combo: str) -> tuple[str, str]:
    parts = combo.split("_")
    if len(parts) != 2 or parts[0] not in MSA_MODES or parts[1] not in TEMPLATE_MODES:
        raise ValueError(
            f"Invalid combination '{combo}' — expected <msaMode>_<templateMode> "
            f"with msaMode in {sorted(MSA_MODES)} and templateMode in {sorted(TEMPLATE_MODES)}."
        )
    return parts[0], parts[1]


def resolve_source_folder(input_dir: Path, msa_mode: str, template_mode: str) -> tuple[Path, set[str]]:
    """Finds the data-generation output folder to build a <msaMode>_<templateMode>
    combination from, and the MSA field(s) that still need to be nulled out of it.

    Data-generation folders are named "msa_<unpaired|paired|full>_<template_onquery|template_standard>".
    An exact <msa_mode> match is preferred — that data was generated with exactly the requested
    MSA type (possibly because the other type was never searched for at all), so nothing needs
    blanking. If none exists, falls back to a "msa_full_" folder and nulls the unwanted field(s),
    which lets several MSA/template combinations be reconstructed from one data-generation run
    without rerunning alignments for each of them.
    """
    template_candidates = (
        list(TEMPLATE_INPUT_SUBFOLDER.values())
        if template_mode == "no"
        else [TEMPLATE_INPUT_SUBFOLDER[template_mode]]
    )

    if msa_mode != "no":
        for template_folder in template_candidates:
            folder = input_dir / f"msa_{msa_mode}_{template_folder}"
            if folder.is_dir():
                return folder, set()

    for template_folder in template_candidates:
        folder = input_dir / f"msa_full_{template_folder}"
        if folder.is_dir():
            return folder, MSA_NULL_FIELDS[msa_mode]

    if msa_mode == "no":
        # Both MSA fields get nulled regardless, so any available MSA type works.
        for template_folder in template_candidates:
            for msa_type in DATAGEN_MSA_TYPES:
                folder = input_dir / f"msa_{msa_type}_{template_folder}"
                if folder.is_dir():
                    return folder, MSA_NULL_FIELDS["no"]

    looked_for = [f"msa_{m}_{t}" for t in template_candidates for m in DATAGEN_MSA_TYPES]
    raise FileNotFoundError(
        f"No data-generation output folder found under {input_dir} for msaMode='{msa_mode}', "
        f"templateMode='{template_mode}' (looked for: {', '.join(looked_for)})"
    )


def write_af_json(base_data, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "alphafold_input.json", "w") as wf:
        json.dump(base_data, wf, indent=4)


def reconstruct_protein(mhc, tra, trb, row, output_path, folder_name, description, null_fields, null_templates):
    id_list = ["D", "E", "C", "A"]
    jsonf = {"name": row["name"]}
    proteins = []
    for id in id_list:
        if id == "D":  # This is the TRA
            sequences = tra.get("sequences")
            if len(sequences) != 1:
                raise ValueError(f"{id}: expected exactly 1 sequence, found {len(sequences)}")
            protein_info = sequences[0].get("protein")
        elif id == "E":  # This is the TRB
            sequences = trb.get("sequences")
            if len(sequences) != 1:
                raise ValueError(f"{id}: expected exactly 1 sequence, found {len(sequences)}")
            protein_info = sequences[0].get("protein")
        elif id == "A":  # This is mhc
            sequences = mhc.get("sequences")
            if len(sequences) != 1:
                raise ValueError(f"{id}: expected exactly 1 sequence, found {len(sequences)}")
            protein_info = sequences[0].get("protein")
        elif id == "C":  # This is peptide — never had MSA/templates to begin with
            proteins.append({
                "protein": {
                    "id": id,
                    "sequence": row["peptide"],
                    "unpairedMsa": "",
                    "pairedMsa": "",
                    "templates": [],
                    "description": description,
                }
            })
            continue

        for field in null_fields:
            protein_info[field] = ""
        if null_templates:
            protein_info["templates"] = []
        proteins.append({"protein": protein_info})

    jsonf["sequences"] = proteins
    jsonf["modelSeeds"] = [1]
    jsonf["dialect"] = "alphafold3"
    jsonf["version"] = 1

    write_af_json(jsonf, output_path / folder_name / row["name"])


def process_protein(args):
    protein_name, folder, row, output_path, folder_name, null_fields, null_templates = args
    mhc_name = row["MHCA_aa_id"]
    tra_name = row["TRA_aa_id"]
    trb_name = row["TRB_aa_id"]
    description = f"{mhc_name}{tra_name}{trb_name}"

    mhc_file = folder / mhc_name / f"{mhc_name}_data.json"
    tra_file = folder / tra_name / f"{tra_name}_data.json"
    trb_file = folder / trb_name / f"{trb_name}_data.json"

    if not (mhc_file.exists() and tra_file.exists() and trb_file.exists()):
        return f"Skipped {protein_name}"

    with open(mhc_file) as f:
        mhc = json.load(f)
    with open(tra_file) as f:
        tra = json.load(f)
    with open(trb_file) as f:
        trb = json.load(f)

    reconstruct_protein(mhc, tra, trb, row, output_path, folder_name, description, null_fields, null_templates)
    return None


def main():
    args = parse_args()
    input_dir = Path(args.input)
    og_data = pd.read_csv(args.data)
    proteins = og_data["name"]

    combos = args.combinations.split()
    parsed_combos = [(combo, *parse_combination(combo)) for combo in combos]
    print(f"Computing reconstruction for {len(parsed_combos)} combinations: {', '.join(combo for combo, _, _ in parsed_combos)}")

    for combo, msa_mode, template_mode in parsed_combos:
        folder, null_fields = resolve_source_folder(input_dir, msa_mode, template_mode)
        null_templates = template_mode == "no"
        folder_name = f"json_{msa_mode}MSA_{template_mode}Template"

        print(f"Reconstructing combination '{msa_mode}MSA_{template_mode}Template' from {folder} -> {args.output / folder_name}")

        tasks = [
            (p, folder, og_data[og_data["name"] == p].iloc[0], args.output, folder_name, null_fields, null_templates)
            for p in proteins
        ]
        with ProcessPoolExecutor(max_workers=10) as ex:
            futures = [ex.submit(process_protein, task) for task in tasks]
            for future in as_completed(futures):
                try:
                    result = future.result()
                    if result is not None:
                        print(result)
                except Exception:
                    print("Worker failed with exception:")
                    raise  # re-raise to crash loudly


if __name__ == "__main__":
    main()
