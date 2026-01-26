import pandas as pd
import argparse
from pathlib import Path
import json
import re

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the csv input file with query data",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output file with the results",
    )
    parser.add_argument(
        "-p",
        "--partition",
        type=int,
        help="Partition for which we want to generate json files, in case we only want specific files",
    )
    return parser.parse_args()


def samplename_to_array(df, output_path, partition = None):
    if partition is not None:
        df = df[df.partition == partition]
    if "name" not in df.columns:
        raise ValueError(
            f"Input file must contain a 'name' column. "
            f"Found columns: {list(df.columns)}"
        )
    df['name'] = df['name'].str.lower()
    sample_names = df['name'].unique()
    output_df = pd.DataFrame({
        'ArrayTaskID': range(1, len(sample_names)+1),
        'SampleName': sample_names
    })

    output_df.to_csv(f"{output_path}/samplename_to_array.txt", sep='\t', index=False)

def get_mhc_sequence(
    allele: str,
) -> str:
    """Returns the amino acid sequence of an MHC allele.

    Args:
        allele: MHC allele name.

    Returns:
        Amino acid sequence of the MHC allele.
    """

    # grab the gene name
    data_dir = Path("databases/mhc_sequences")
    gene = re.search(r"([A-Z])\*.*", allele).group(1)
    with open(data_dir / f"{gene}_prot.json") as f:
        data = json.load(f)

    # default to 01 for missing fields
    num_fields = len(allele.split(":"))
    allele = allele + ":01" * (4 - num_fields)

    # grab the sequence
    try:
        seq = data[allele]
    except KeyError:
        raise KeyError(f"Allele {allele} not found in {data_dir}")
    return seq

def json_af3(df, output, partition = None):

    if partition is not None:
        df = df[df.partition == partition]

    df['name'] = df['name'].str.lower()
    
    if "epitope_aa" in df.columns:
        df = df.rename(columns={"epitope_aa": "peptide"})
        
    id_list = ["D", "E", "C", "A"]
    for idx, row in df.iterrows():
        jsonf = {}
        jsonf["name"] = row["name"]

        seq_list = [row["TRA_aa"], row["TRB_aa"], row["peptide"], row["MHCA_aa"]]
        proteins = []
        for i in range(len(id_list)):
            if id_list[i] == "C":
                protein = {
                    "protein": {
                        "id": id_list[i],
                        "sequence": seq_list[i],
                        "unpairedMsa": "",
                        "pairedMsa": "",
                        "templates": []

                    }
                }

                proteins.append(protein)
            else: 
                protein = {
                    "protein": {
                        "id": id_list[i],
                        "sequence": seq_list[i]
                    }
                }
                
                proteins.append(protein)

        jsonf["sequences"] = proteins
        jsonf["modelSeeds"] = [1]
        jsonf["dialect"] = "alphafold3"
        jsonf["version"] = 1

        
        output_dir = output / "json_msa_template" / row["name"]
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(output_dir / "alphafold_input.json", "w") as f:
            json.dump(jsonf, f, indent=4)


def main() -> None:
    """Main entry point to modeling pipeline."""
    args = parse_args()
    input_csv = args.input
    output_path = args.output
    output_path.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_csv)
    if "name" not in df.columns:
        df["name"] = df["peptide"] + "_" + df["A1"] + "_" + df["A2"] + "_" + df["A3"] + "_" + df["B1"] + "_" + df["B2"] + "_" + df["B3"] 
    samplename_to_array(df, output_path, partition = args.partition)
    if "MHCA_aa" not in df.columns:
        df["MHCA_aa"] = df["allele"].map(
            {allele: get_mhc_sequence(allele) for allele in df["allele"].unique()}
        )
        df["MHCB_aa"] = ""
    df.to_csv(Path(output_path/ f"{Path(input_csv).stem}_hla.csv"))
    json_af3(df, Path(output_path / "jsonFiles/"), partition = args.partition)

if __name__ == "__main__":
    main()


