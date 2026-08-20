#Validates create_custom_json.py output against its own source data: for every
#json_<msaMode>MSA_<templateMode>Template/<name>/alphafold_input.json, checks that
#pairedMsa/unpairedMsa on the TRA/TRB/MHC chains are blanked exactly where the folder's
#msaMode says they should be (and left untouched otherwise), and that the templates field
#matches the exact source data-generation folder used (or is emptied for templateMode=no).
#Usage example: python validate_custom_json.py -i .../dataPipelineOut -o .../customJSON -d data_hla_withid.csv
import argparse
import json
import re
from pathlib import Path
import pandas as pd

from structureTCR.jsonPrep.create_custom_json import (
    MSA_NULL_FIELDS,
    TEMPLATE_INPUT_SUBFOLDER,
    resolve_template_folder,
)

FOLDER_NAME_RE = re.compile(r"^json_(?P<msa>unpaired|paired|full|no)MSA_(?P<template>onquery|standard|no)Template$")
MSA_FIELDS = {"pairedMsa", "unpairedMsa"}
CHAIN_ID_COLUMN = {"D": "TRA_aa_id", "E": "TRB_aa_id", "A": "MHCA_aa_id"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validates create_custom_json.py output against its source data-generation folders."
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True,
        help="Path to the data-generation output (parent of template_onquery/template_standard) — same as create_custom_json.py's -i",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True,
        help="Path to the customJSON folder to validate (parent of json_<combination> folders) — same as create_custom_json.py's -o",
    )
    parser.add_argument(
        "-d", "--data", type=Path, required=True,
        help="Path to csv with original datapoints and unique sequence identifiers — same as create_custom_json.py's -d",
    )
    return parser.parse_args()


def load_chain_protein(source_root: Path, template_mode: str, chain_id: str):
    folder = resolve_template_folder(source_root, template_mode)
    path = folder / chain_id / f"{chain_id}_data.json"
    if not path.exists():
        return folder, None
    with open(path) as f:
        data = json.load(f)
    return folder, data["sequences"][0]["protein"]


def check_protein(json_path: Path, row, msa_mode: str, template_mode: str, source_root: Path) -> list[str]:
    with open(json_path) as f:
        data = json.load(f)

    null_fields = MSA_NULL_FIELDS[msa_mode]
    kept_fields = MSA_FIELDS - null_fields
    errors = []

    for seq in data["sequences"]:
        prot = seq["protein"]
        cid = prot["id"]
        if cid == "C":  # peptide has no source file, its fields are always fixed
            continue

        chain_id = row[CHAIN_ID_COLUMN[cid]]
        folder, source_prot = load_chain_protein(source_root, template_mode, chain_id)
        if source_prot is None:
            errors.append(f"{cid}: no source data found for chain {chain_id} under {folder}")
            continue

        for field in null_fields:
            if prot.get(field, None) != "":
                errors.append(f"{cid}: expected {field} to be blanked (msaMode={msa_mode}), got {len(prot.get(field) or '')} chars")
        for field in kept_fields:
            if prot.get(field) != source_prot.get(field):
                errors.append(f"{cid}: {field} does not match source (msaMode={msa_mode} should leave it untouched)")

        if template_mode == "no":
            if prot.get("templates"):
                errors.append(f"{cid}: expected templates=[] (templateMode=no), found {len(prot['templates'])}")
        else:
            if prot.get("templates") != source_prot.get("templates"):
                errors.append(
                    f"{cid}: templates do not match source {TEMPLATE_INPUT_SUBFOLDER[template_mode]}/ "
                    f"(output has {len(prot.get('templates') or [])}, source has {len(source_prot.get('templates') or [])})"
                )

    return errors


def main():
    args = parse_args()
    og_data = pd.read_csv(args.data)

    combo_folders = sorted(p for p in args.output.iterdir() if p.is_dir() and FOLDER_NAME_RE.match(p.name))
    if not combo_folders:
        raise SystemExit(f"No json_<msaMode>MSA_<templateMode>Template folders found under {args.output}")

    total_errors = 0
    total_checked = 0
    for folder in combo_folders:
        m = FOLDER_NAME_RE.match(folder.name)
        msa_mode, template_mode = m["msa"], m["template"]
        protein_dirs = sorted(p for p in folder.iterdir() if p.is_dir())
        print(f"Checking {folder.name} ({len(protein_dirs)} proteins)...")

        for protein_dir in protein_dirs:
            json_path = protein_dir / "alphafold_input.json"
            if not json_path.exists():
                print(f"  {protein_dir.name}: MISSING alphafold_input.json")
                total_errors += 1
                continue

            rows = og_data[og_data["name"] == protein_dir.name]
            if rows.empty:
                print(f"  {protein_dir.name}: not found in {args.data}, skipping")
                continue

            total_checked += 1
            errors = check_protein(json_path, rows.iloc[0], msa_mode, template_mode, args.input)
            for e in errors:
                print(f"  {protein_dir.name}: {e}")
                total_errors += 1

    print()
    if total_errors:
        print(f"FAILED: {total_errors} issue(s) found across {total_checked} protein(s) checked")
        raise SystemExit(1)
    print(f"All checks passed ({total_checked} protein(s) across {len(combo_folders)} folder(s)).")


if __name__ == "__main__":
    main()
