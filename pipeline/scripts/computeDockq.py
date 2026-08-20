from Bio.PDB import PDBIO, MMCIFParser
import os
import subprocess
from pathlib import Path
from typing import Union
import argparse
import json
from concurrent.futures import ProcessPoolExecutor


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the folder with af3 output for which to compute the DockQ score",
    )
    parser.add_argument(
        "-t",
        "--templates",
        type=Path,
        required=True,
        help="Path to the folder with model templates in pdb format",
    )
    parser.add_argument(
        "-n1",
        "--native_chain1",
        nargs='+',
        #required = True,
        help="Chain ID of first native chain",
    )
    parser.add_argument(
        "-n2",
        "--native_chain2",
        nargs='+',
        #required = True,
        help="Chain ID of second native chain",
    )
    parser.add_argument(
        "-m1",
        "--model_chain1",
        nargs='+',
        #required = True,
        help="Chain ID of first model chain",
    )
    parser.add_argument(
        "-m2",
        "--model_chain2",
        nargs='+',
        #required = True,
        help="Chain ID of second model chain",
    )
    parser.add_argument(
        "-r", 
        "--recompute", 
        nargs='+',
        help="List of proteins for which to recompute the metrics, if not specified it will compute the metrics for all the complexes",
        )
    parser.add_argument(
        "-s", 
        "--suffix", 
        help="Specify molecules used to compute dockq",
        )

    return parser.parse_args()

def cif_to_pdb(path_f, pdb_id, model_name):
    #print("Converting cif to pdb")
    cif_parser = MMCIFParser(QUIET=True)
    sample = path_f /f"{pdb_id}_{model_name}_model.cif"
    model = cif_parser.get_structure(
        "-",
        open(sample, "rt")
    )
    
    p = PDBIO()
    p.set_structure(model)
    p.save(f"{os.path.splitext(sample)[0]}.pdb")
    return f"{os.path.splitext(sample)[0]}.pdb"


def compute_dockq(predicted_path: Path, 
                  native_path:Path, 
                  native_chain1:  Union[str, list], 
                  model_chain1: Union[str, list], 
                  native_chain2: Union[str, list], 
                  model_chain2: Union[str, list]):
    
    """Runs DockQ on a predicted and ground truth structure.

    Args:
        dockq_path: Path to DockQ.py script.
        predicted_path: Path to predicted structure.
        native_path: Path to ground truth structure.
        native_chain1: Chain ID of first native chain.
        native_chain2: Chain ID of second native chain.
        model_chain1: Chain ID of first model chain.
        model_chain2: Chain ID of second model chain.

    Returns:
        DockQ output.
    """
    #print("Computing DockQ")
    if not isinstance(native_chain1, str):
        native_chain1 = " ".join(native_chain1)
    if not isinstance(native_chain2, str):
        native_chain2 = " ".join(native_chain2)
    if not isinstance(model_chain1, str):
        model_chain1 = " ".join(model_chain1)
    if not isinstance(model_chain2, str):
        model_chain2 = " ".join(model_chain2)
    dockqPath="/mnt/dockq_repo/DockQ.py"
    cmd = f"python3 {dockqPath} {str(predicted_path)} {str(native_path)} -native_chain1 {native_chain1} -model_chain1 {model_chain1} -native_chain2 {native_chain2} -model_chain2 {model_chain2} -no_needle"
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, encoding="utf-8")

def parse_dockq_output(dockq_output: str, output_json_path:Path, pdb_id, model_name, suffix) -> dict:
    """Parses the output of DockQ.

    Args:
        dockq_output: Output string from DockQ.

    Returns:
        Dictionary with DockQ score, fnat, iRMS, and LRMS.
    """
    lines = dockq_output.strip().split("\n")
    results = {}
    for line in lines:
        if line.startswith("DockQ"):
            results["DockQ"] = float(line.split()[-1])
        elif line.startswith("Fnat"):
            results["Fnat"] = float(line.split()[1])
        elif line.startswith("Fnonnat"):
            results["Fnonnat"] = float(line.split()[1])
        elif line.startswith("iRMS"):
            results["iRMS"] = float(line.split()[-1])
        elif line.startswith("LRMS"):
            results["LRMS"] = float(line.split()[-1])
    
    with open(f"{output_json_path}/dockQ_metrics_{pdb_id}_{model_name}{suffix}.json", 'w') as f:
        json.dump(results, f, indent=4)

    return results


def dockq_one_complex(predicted_path, native_path, native_chain1, model_chain1, native_chain2, model_chain2, pdb_id, model_name, suffix):
    pdb_path = cif_to_pdb(predicted_path, pdb_id, model_name)
    result = compute_dockq(pdb_path, native_path, native_chain1, model_chain1, native_chain2, model_chain2)
    output = parse_dockq_output(result.stdout, os.path.splitext(predicted_path)[0], pdb_id, model_name, suffix)

def process_complex(args):
    c, args_ns, suffix = args
    pdb_id = c.name

    template_path = args_ns.templates / f"{pdb_id}.trunc.fit.pdb"
    models = [m for m in c.iterdir() if m.is_dir()]
    for m in models:
        model_name = m.name
        dockq_one_complex(
            predicted_path=m,
            native_path=template_path,
            model_chain1=args_ns.model_chain1,
            native_chain1=args_ns.native_chain1,
            model_chain2=args_ns.model_chain2,
            native_chain2=args_ns.native_chain2,
            pdb_id=pdb_id,
            model_name=model_name,
            suffix = suffix
        )

    return pdb_id

def main():
    args = parse_args()
    input_dir = Path(args.input)
    complexes = [f for f in input_dir.iterdir() if f.is_dir()]
    suffix = args.suffix
    if args.recompute:
        complexes = [c for c in complexes if c.name in args.recompute]

    tasks = [(c, args, suffix) for c in complexes]

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(process_complex, tasks)


if __name__ == "__main__":
    main()