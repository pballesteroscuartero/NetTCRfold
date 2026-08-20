#Extracts the PDB template IDs used in customJSON output for one or more folders
#(json_<combination>/<pdb_id>/alphafold_input.json), producing a long-format dataframe:
#pdb_id, template, folder — one row per template entry found, across all chains (TRA/TRB/MHC).
#Useful for comparing which templates got selected between two combinations
#(e.g. onquery vs standard template mode) for the same datapoints.
#Usage example: python compare_templates.py -f .../json_fullMSA_onqueryTemplate .../json_fullMSA_standardTemplate -o templates.csv
import argparse
import json
import re
from pathlib import Path
import pandas as pd

MMCIF_ID_RE = re.compile(r"(data_\S+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extracts template PDB IDs used in customJSON output, for comparing template selection across folders."
    )
    parser.add_argument(
        "-f",
        "--folders",
        type=Path,
        nargs="+",
        required=True,
        help="Two or more json_<combination> folders (each containing <pdb_id>/alphafold_input.json) to extract templates from",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Optional path to write the resulting dataframe as csv",
    )
    return parser.parse_args()


def extract_mmcif_id(mmcif: str) -> str:
    match = MMCIF_ID_RE.match(mmcif)
    if not match:
        raise ValueError(f"Could not find a 'data_XXX' identifier at the start of mmcif block: {mmcif[:50]!r}")
    return match.group(1)


def extract_templates(folder: Path) -> pd.DataFrame:
    rows = []
    for protein_dir in sorted(p for p in folder.iterdir() if p.is_dir()):
        json_path = protein_dir / "alphafold_input.json"
        if not json_path.exists():
            continue
        with open(json_path) as f:
            data = json.load(f)
        for seq in data["sequences"]:
            prot = seq["protein"]
            for template in prot.get("templates") or []:
                rows.append({
                    "pdb_id": protein_dir.name,
                    "template": extract_mmcif_id(template["mmcif"]),
                    "folder": folder.name,
                })
    return pd.DataFrame(rows, columns=["pdb_id", "template", "folder"])


def main():
    args = parse_args()
    dfs = [extract_templates(folder) for folder in args.folders]
    df = pd.concat(dfs, ignore_index=True)

    print(df)
    print(f"\n{len(df)} template rows across {len(args.folders)} folder(s)")
    for folder_name, sub in df.groupby("folder"):
        print(f"  {folder_name}: {sub['pdb_id'].nunique()} pdb_id(s), {len(sub)} template rows")

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
