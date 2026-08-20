# ipsae.py
# script for calculating the ipSAE score for scoring pairwise protein-protein interactions in AlphaFold2 and AlphaFold3 models
# https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2

# Also calculates:
#    pDockQ: Bryant, Pozotti, and Eloffson. https://www.nature.com/articles/s41467-022-28865-w
#    pDockQ2: Zhu, Shenoy, Kundrotas, Elofsson. https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714
#    LIS: Kim, Hu, Comjean, Rodiger, Mohr, Perrimon. https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1

# Roland Dunbrack
# Fox Chase Cancer Center
# version 4
# January 3, 2026: Fixed Boltz2 issues (PDB and mmCIF format; chainIDs)
# MIT license: script can be modified and redistributed for non-commercial and commercial use, as long as this information is reproduced.

# includes support for Boltz structures and structures with nucleic acids

# It may be necessary to install numpy with the following command:
#      pip install numpy

# Usage:

#  python ipsae.py <path_to_af2_pae_file>        <path_to_af2_pdb_file>     <pae_cutoff> <dist_cutoff>
#  python ipsae.py <path_to_af3_pae_file>        <path_to_af3_cif_file>     <pae_cutoff> <dist_cutoff>
#  python ipsae.py <path_to_boltz_pae_npz_file>  <path_to_boltz_cif_file>   <pae_cutoff> <dist_cutoff>
#  python ipsae.py <path_to_boltz_pae_npz_file>  <path_to_boltz_pdb_file>   <pae_cutoff> <dist_cutoff>
#
# All output files will be in same path/folder as cif or pdb file

import sys, os, math
import json
import numpy as np
import pandas as pd

np.set_printoptions(threshold=np.inf)  # for printing out full numpy arrays for debugging

def ptm_func(x,d0):
    return 1.0/(1+(x/d0)**2.0)

def calc_d0(L,pair_type):
    L=float(L)
    min_value=1.0
    if pair_type=='nucleic_acid': min_value=2.0
    if L>27:
        d0=1.24*(L-15)**(1.0/3.0) - 1.8
    else:
        d0=1.0
    return max(min_value, d0)


# Define the d0 functions for numbers and arrays; minimum value = 1.0; from Yang and Skolnick, PROTEINS: Structure, Function, and Bioinformatics 57:702–710 (2004)
def calc_d0_array(L, pair_type):
    # Convert L to a NumPy array if it isn't already one (enables flexibility in input types)
    # fixed 01.03.2026: now returns 1.00 instead of 1.04 for minimum value
    L = np.array(L, dtype=float)
    L = np.maximum(26,L)
    min_value=1.0

    if pair_type=='nucleic_acid': min_value=2.0

    # Calculate d0 using the vectorized operation
    return np.maximum(min_value, 1.24 * (L - 15) ** (1.0/3.0) - 1.8)

# Define the parse_atom_line function for PDB lines (by column) and mmCIF lines (split by white_space)
# parsed_line = parse_atom_line(line)
# line = "ATOM    123  CA  ALA A  15     11.111  22.222  33.333  1.00 20.00           C"
def parse_pdb_atom_line(line):
    atom_num = line[6:11].strip()
    atom_name = line[12:16].strip()
    residue_name = line[17:20].strip()
    if residue_name == "LIG": return None  # ligands in Boltz PDB-format files
    chain_id = line[21].strip()
    residue_seq_num = line[22:26].strip()
    x = line[30:38].strip()
    y = line[38:46].strip()
    z = line[46:54].strip()

    # Convert string numbers to integers or floats as appropriate
    atom_num = int(atom_num)
    residue_seq_num = int(residue_seq_num)
    x = float(x)
    y = float(y)
    z = float(z)

    return {
        'atom_num': atom_num,
        'atom_name': atom_name,
        'residue_name': residue_name,
        'chain_id': chain_id,
        'residue_seq_num': residue_seq_num,
        'x': x,
        'y': y,
        'z': z
    }

def parse_cif_atom_line(line,fielddict):
    # for parsing AF3 and Boltz1/2 mmCIF files
    # ligands do not have residue numbers but modified residues do. Return "None" for ligand.
    # AF3 mmcif lines
    # 0      1   2   3     4  5  6 7  8  9  10      11     12      13   14    15 16 17
    #ATOM   1294 N   N     . ARG A 1 159 ? 5.141   -14.096 10.526  1.00 95.62 159 A 1
    #ATOM   1295 C   CA    . ARG A 1 159 ? 4.186   -13.376 11.366  1.00 96.27 159 A 1
    #ATOM   1296 C   C     . ARG A 1 159 ? 2.976   -14.235 11.697  1.00 96.42 159 A 1
    #ATOM   1297 O   O     . ARG A 1 159 ? 2.654   -15.174 10.969  1.00 95.46 159 A 1
    # ...
    #HETATM 1305 N   N     . TPO A 1 160 ? 2.328   -13.853 12.742  1.00 96.42 160 A 1
    #HETATM 1306 C   CA    . TPO A 1 160 ? 1.081   -14.560 13.218  1.00 96.78 160 A 1
    #HETATM 1307 C   C     . TPO A 1 160 ? -2.115  -11.668 12.263  1.00 96.19 160 A 1
    #HETATM 1308 O   O     . TPO A 1 160 ? -1.790  -11.556 11.113  1.00 95.75 160 A 1
    # ...
    #HETATM 2608 P   PG    . ATP C 3 .   ? -6.858  4.182   10.275  1.00 84.94 1   C 1
    #HETATM 2609 O   O1G   . ATP C 3 .   ? -6.178  5.238   11.074  1.00 75.56 1   C 1
    #HETATM 2610 O   O2G   . ATP C 3 .   ? -5.889  3.166   9.748   1.00 75.15 1   C 1
    # ...
    #HETATM 2639 MG  MG    . MG  D 4 .   ? -7.262  2.709   4.825   1.00 91.47 1   D 1
    #HETATM 2640 MG  MG    . MG  E 5 .   ? -4.994  2.251   8.755   1.00 85.96 1   E 1


    # Boltz1 mmcif files (in non-standard order))
    #_atom_site.group_PDB
    #_atom_site.id
    #_atom_site.type_symbol
    #_atom_site.label_atom_id
    #_atom_site.label_alt_id
    #_atom_site.label_comp_id
    #_atom_site.label_seq_id
    #_atom_site.auth_seq_id
    #_atom_site.pdbx_PDB_ins_code
    #_atom_site.label_asym_id
    #_atom_site.Cartn_x
    #_atom_site.Cartn_y
    #_atom_site.Cartn_z
    #_atom_site.occupancy
    #_atom_site.label_entity_id
    #_atom_site.auth_asym_id
    #_atom_site.auth_comp_id
    #_atom_site.B_iso_or_equiv
    #_atom_site.pdbx_PDB_model_num
    # 0     1     2  3     4   5   6    7   8  9  10          11         12       13  14 15 16  17 18
    #ATOM   2652  N  N     . ASN  43   43   ?  B  10.83538   6.06359    18.45139   1  2  B  ASN  1  1
    #ATOM   2653  C  CA    . ASN  43   43   ?  B  10.76295   5.07366    19.53232   1  2  B  ASN  1  1
    #ATOM   2654  C  C     . ASN  43   43   ?  B  11.21770   5.64437    20.88774   1  2  B  ASN  1  1
    #ATOM   2655  O  O     . ASN  43   43   ?  B  12.06730   6.51688    20.91168   1  2  B  ASN  1  1
    #ATOM   2656  C  CB    . ASN  43   43   ?  B  11.60137   3.84778    19.19481   1  2  B  ASN  1  1
    #ATOM   2657  C  CG    . ASN  43   43   ?  B  10.96208   3.03997    18.07013   1  2  B  ASN  1  1
    #ATOM   2658  O  OD1   . ASN  43   43   ?  B  9.79094    3.17033    17.81165   1  2  B  ASN  1  1
    #ATOM   2659  N  ND2   . ASN  43   43   ?  B  11.77101   2.23791    17.39764   1  2  B  ASN  1  1
    #HETATM 2660  P  PG    . ATP  .    1    ?  C  -8.79525   6.04621    -4.99212   1  3  C  ATP  1  1
    #HETATM 2661  O  O1G   . ATP  .    1    ?  C  -10.01901  6.83468    -5.24825   1  3  C  ATP  1  1
    #HETATM 2662  O  O2G   . ATP  .    1    ?  C  -9.03047   4.56941    -4.85246   1  3  C  ATP  1  1
    #HETATM 2663  O  O3G   . ATP  .    1    ?  C  -7.97335   6.60305    -3.86656   1  3  C  ATP  1  1
    #HETATM 2664  P  PB    . ATP  .    1    ?  C  -6.63618   7.04315    -6.56073   1  3  C  ATP  1  1
    #HETATM 2665  O  O1B   . ATP  .    1    ?  C  -7.04640   8.36577    -7.14326   1  3  C  ATP  1  1
    #HETATM 2666  O  O2B   . ATP  .    1    ?  C  -5.79036   7.13926    -5.33995   1  3  C  ATP  1  1


    linelist =        line.split()
    atom_num =        linelist[ fielddict['id'] ]
    atom_name =       linelist[ fielddict['label_atom_id'] ]
    residue_name =    linelist[ fielddict['label_comp_id'] ]
    if "auth_asym_id" in fielddict:
        chain_id =    linelist[ fielddict['auth_asym_id'] ]
    else:
        chain_id =    linelist[ fielddict['label_asym_id'] ]

    residue_seq_num = linelist[ fielddict['label_seq_id'] ]
    x =               linelist[ fielddict['Cartn_x'] ]
    y =               linelist[ fielddict['Cartn_y'] ]
    z =               linelist[ fielddict['Cartn_z'] ]

    if residue_seq_num == ".": return None   # ligand atom

    # Convert string numbers to integers or floats as appropriate
    atom_num = int(atom_num)
    residue_seq_num = int(residue_seq_num)
    x = float(x)
    y = float(y)
    z = float(z)

    return {
        'atom_num': atom_num,
        'atom_name': atom_name,
        'residue_name': residue_name,
        'chain_id': chain_id,
        'residue_seq_num': residue_seq_num,
        'x': x,
        'y': y,
        'z': z
    }


# Function for printing out residue numbers in PyMOL scripts
def contiguous_ranges(numbers):
    if not numbers:  # Check if the set is empty
        return

    sorted_numbers = sorted(numbers)  # Sort the numbers
    start = sorted_numbers[0]
    end = start
    ranges = []  # List to store ranges

    def format_range(start, end):
        if start == end:
            return f"{start}"
        else:
            return f"{start}-{end}"

    for number in sorted_numbers[1:]:
        if number == end + 1:
            end = number
        else:
            ranges.append(format_range(start, end))
            start = end = number

    # Append the last range after the loop
    ranges.append(format_range(start, end))

    # Join all ranges with a plus sign and print the result
    string='+'.join(ranges)
    return(string)

# Initializes a nested dictionary with all values set to 0
def init_chainpairdict_zeros(chainlist):
    return {chain1: {chain2: 0 for chain2 in chainlist if chain1 != chain2} for chain1 in chainlist}

# Initializes a nested dictionary with NumPy arrays of zeros of a specified size
def init_chainpairdict_npzeros(chainlist, arraysize):
    return {chain1: {chain2: np.zeros(arraysize) for chain2 in chainlist if chain1 != chain2} for chain1 in chainlist}

# Initializes a nested dictionary with empty sets.
def init_chainpairdict_set(chainlist):
    return {chain1: {chain2: set() for chain2 in chainlist if chain1 != chain2} for chain1 in chainlist}

def classify_chains(chains, residue_types):
    nuc_residue_set = {"DA", "DC", "DT", "DG", "A", "C", "U", "G"}
    chain_types = {}

    # Get unique chains and iterate over them
    _, first_idx = np.unique(chains, return_index=True)
    unique_chains = chains[np.sort(first_idx)]

    for chain in unique_chains:
        # Find indices where the current chain is located
        indices = np.where(chains == chain)[0]
        # Get the residues for these indices
        chain_residues = residue_types[indices]
        # Count nucleic acid residues
        nuc_count = sum(residue in nuc_residue_set for residue in chain_residues)

        # Determine if the chain is a nucleic acid or protein
        chain_types[chain] = 'nucleic_acid' if nuc_count > 0 else 'protein'

    return chain_types

import numpy as np

def build_ipsae_inputs_from_structure(structure, residue_set):
    """
    Replicates ipsae.py's manual CIF/PDB line-parsing block, but builds
    residues / cb_residues / chains / token_mask directly from a biotite
    AtomArray that's already been loaded (e.g. via strucio.load_structure),
    instead of re-opening and parsing the file line-by-line.

    Requires: structure.atom_id must already be set (e.g. via
    structure.set_annotation("atom_id", np.arange(1, len(structure)+1))
    as you already do in process_complex).

    IMPORTANT: verify the `is_ligand` logic below against your actual data
    before trusting this in production -- see note under the loop.
    """
    residues = []
    cb_residues = []
    chains = []
    token_mask = []

    atom_names = structure.atom_name
    res_names  = structure.res_name
    chain_ids  = structure.chain_id
    res_ids    = structure.res_id
    atom_ids   = structure.atom_id
    coords     = structure.coord
    hetero     = structure.hetero if hasattr(structure, "hetero") else np.zeros(len(structure), dtype=bool)

    for i in range(len(structure)):
        atom_name = atom_names[i]
        res_name  = res_names[i]
        chain_id  = chain_ids[i]
        res_id    = res_ids[i]

        # Original logic: CIF ligands have auth_seq_id == "." (no residue number);
        # PDB ligands have residue_name == "LIG". Biotite doesn't preserve the
        # "." literal, so this approximates it via hetero flag + missing/invalid
        # residue id. VALIDATE this against a known ligand-containing structure
        # before relying on it -- it's the one part of this rewrite that isn't
        # a direct 1:1 translation.
        is_ligand = hetero[i] and (res_id is None or res_id <= 0)

        if is_ligand:
            token_mask.append(0)
            continue

        if atom_name == "CA" or "C1" in atom_name:
            token_mask.append(1)
            residues.append({
                'atom_num': int(atom_ids[i]),
                'coor': coords[i].copy(),
                'res': res_name,
                'chainid': chain_id,
                'resnum': int(res_id),
                'residue': f"{res_name:3}   {chain_id:3} {res_id:4}"
            })
            chains.append(chain_id)

        if atom_name == "CB" or "C3" in atom_name or (res_name == "GLY" and atom_name == "CA"):
            cb_residues.append({
                'atom_num': int(atom_ids[i]),
                'coor': coords[i].copy(),
                'res': res_name,
                'chainid': chain_id,
                'resnum': int(res_id),
                'residue': f"{res_name:3}   {chain_id:3} {res_id:4}"
            })

        if atom_name != "CA" and "C1" not in atom_name and res_name not in residue_set:
            token_mask.append(0)

    return residues, cb_residues, chains, token_mask


def compute_ipsae(pae_data, pae_data_summary, pdb_file, pdb_path, pae_cutoff, dist_cutoff, write_file = False):
    """
    Computes the ipSAE score for a given PAE matrix and PDB file.
    pae_file_path: Path to the PAE file. It's the _confidences.json file
    pdb_path: Path to the PDB file. It's the .pdb or .cif file output from AF3 modeling
    pae_cutoff: The PAE cutoff value for calculating ipSAE.
    dist_cutoff: The distance cutoff value for calculating ipSAE.
    """

    pae_string = str(int(pae_cutoff))
    if pae_cutoff<10:  pae_string="0"+pae_string
    dist_string =      str(int(dist_cutoff))
    if dist_cutoff<10: dist_string="0"+dist_string
    af3 =    True
    #file2_path =       path_stem + "_byres.txt"
    #pml_path =         path_stem + ".pml"

    
    #PML =              open(pml_path,'w')
    #OUT2 =             open(file2_path,'w')

    ptm_func_vec=np.vectorize(ptm_func)  # vector version

    residue_set= {"ALA", "ARG", "ASN", "ASP", "CYS",
              "GLN", "GLU", "GLY", "HIS", "ILE",
              "LEU", "LYS", "MET", "PHE", "PRO",
              "SER", "THR", "TRP", "TYR", "VAL",
              "DA", "DC", "DT", "DG", "A", "C", "U", "G"}
    residues, cb_residues, chains, token_mask = build_ipsae_inputs_from_structure(pdb_file, residue_set)
    numres = len(residues)
    CA_atom_num=  np.array([res['atom_num']-1 for res in residues])  # for AF3 atom indexing from 0
    CB_atom_num=  np.array([res['atom_num']-1 for res in cb_residues])  # for AF3 atom indexing from 0
    coordinates = np.array([res['coor']       for res in cb_residues])
    chains = np.array(chains)

    _, first_idx = np.unique(chains, return_index=True)
    unique_chains = chains[np.sort(first_idx)]
    token_array=np.array(token_mask)
    ntokens=np.sum(token_array)
    residue_types=np.array([res['res'] for res in residues])

    # chain types (nucleic acid (NA) or protein) and chain_pair_types ('nucleic_acid' if either chain is NA) for d0 calculation
    # arbitrarily setting d0 to 2.0 for NA/protein or NA/NA chain pairs (approximately 21 base pairs)
    d0_nucleic_acid=2.0
    chain_dict = classify_chains(chains, residue_types)
    chain_pair_type = init_chainpairdict_zeros(unique_chains)
    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1==chain2: continue
            if chain_dict[chain1] == 'nucleic_acid' or chain_dict[chain2] == 'nucleic_acid':
                chain_pair_type[chain1][chain2]='nucleic_acid'
            else:
                chain_pair_type[chain1][chain2]='protein'

    # Calculate distance matrix using NumPy broadcasting
    distances = np.sqrt(((coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :])**2).sum(axis=2))


    if "atom_plddts" in pae_data:
        atom_plddts=np.array(pae_data['atom_plddts'])
        plddt=atom_plddts[CA_atom_num]  # pull out residue plddts from Calpha atoms
    else:
        plddt = np.zeros(numres)

    # Get pairwise residue PAE matrix by identifying one token per protein residue.
    # Modified residues have separate tokens for each atom, so need to pull out Calpha atom as token
    # Skip ligands
    if 'pae' in pae_data:
        pae_matrix_af3 = np.array(pae_data['pae'])
    else:
        print("no PAE data in AF3 json file; quitting")
        sys.exit()

    # Set pae_matrix for AF3 from subset of full PAE matrix from json file
    token_array=np.array(token_mask)
    pae_matrix = pae_matrix_af3[np.ix_(token_array.astype(bool), token_array.astype(bool))]
    # Get iptm matrix from AF3 summary_confidences file
    iptm_af3=   {chain1: {chain2: 0     for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

    af3_chain_pair_iptm_data=pae_data_summary['chain_pair_iptm']
    for nchain1, chain1 in enumerate(unique_chains):
        for nchain2, chain2 in enumerate(unique_chains):
            if chain1 == chain2: continue
            iptm_af3[chain1][chain2]=af3_chain_pair_iptm_data[nchain1][nchain2]
    
    # Compute chain-pair-specific interchain PTM and PAE, count valid pairs, and count unique residues
    # First, create dictionaries of appropriate size: top keys are chain1 and chain2 where chain1 != chain2
    # Nomenclature:
    # iptm_d0chn =  calculate iptm  from PAEs with no PAE cutoff; d0 = numres in chain pair = len(chain1) + len(chain2)
    # ipsae_d0chn = calculate ipsae from PAEs with PAE cutoff;    d0 = numres in chain pair = len(chain1) + len(chain2)
    # ipsae_d0dom = calculate ipsae from PAEs with PAE cutoff;    d0 from number of residues in chain1 and chain2 that have interchain PAE<cutoff
    # ipsae_d0res = calculate ipsae from PAEs with PAE cutoff;    d0 from number of residues in chain2 that have interchain PAE<cutoff given residue in chain1
    #
    # for each chain_pair iptm/ipsae, there is (for example)
    # ipsae_d0res_byres = by-residue array;
    # ipsae_d0res_asym  = asymmetric pair value (A->B is different from B->A)
    # ipsae_d0res_max   = maximum of A->B and B->A value
    # ipsae_d0res_asymres = identify of residue that provides each asym maximum
    # ipsae_d0res_maxres =  identify of residue that provides each maximum over both chains
    #
    # n0num = number of residues in whole complex provided by AF2 model
    # n0chn = number of residues in chain pair = len(chain1) + len(chain2)
    # n0dom = number of residues in chain pair that have good PAE values (<cutoff)
    # n0res = number of residues in chain2 that have good PAE residues for each residue of chain1

    iptm_d0chn_byres  = init_chainpairdict_npzeros(unique_chains, numres)
    ipsae_d0chn_byres = init_chainpairdict_npzeros(unique_chains, numres)
    ipsae_d0dom_byres = init_chainpairdict_npzeros(unique_chains, numres)
    ipsae_d0res_byres = init_chainpairdict_npzeros(unique_chains, numres)

    iptm_d0chn_asym   = init_chainpairdict_zeros(unique_chains)
    ipsae_d0chn_asym  = init_chainpairdict_zeros(unique_chains)
    ipsae_d0dom_asym  = init_chainpairdict_zeros(unique_chains)
    ipsae_d0res_asym  = init_chainpairdict_zeros(unique_chains)

    iptm_d0chn_max    = init_chainpairdict_zeros(unique_chains)
    ipsae_d0chn_max   = init_chainpairdict_zeros(unique_chains)
    ipsae_d0dom_max   = init_chainpairdict_zeros(unique_chains)
    ipsae_d0res_max   = init_chainpairdict_zeros(unique_chains)

    iptm_d0chn_asymres   = init_chainpairdict_zeros(unique_chains)
    ipsae_d0chn_asymres  = init_chainpairdict_zeros(unique_chains)
    ipsae_d0dom_asymres  = init_chainpairdict_zeros(unique_chains)
    ipsae_d0res_asymres  = init_chainpairdict_zeros(unique_chains)

    iptm_d0chn_maxres    = init_chainpairdict_zeros(unique_chains)
    ipsae_d0chn_maxres   = init_chainpairdict_zeros(unique_chains)
    ipsae_d0dom_maxres   = init_chainpairdict_zeros(unique_chains)
    ipsae_d0res_maxres   = init_chainpairdict_zeros(unique_chains)

    n0chn       = init_chainpairdict_zeros(unique_chains)
    n0dom       = init_chainpairdict_zeros(unique_chains)
    n0dom_max   = init_chainpairdict_zeros(unique_chains)
    n0res       = init_chainpairdict_zeros(unique_chains)
    n0res_max   = init_chainpairdict_zeros(unique_chains)
    n0res_byres = init_chainpairdict_npzeros(unique_chains, numres)

    d0chn       = init_chainpairdict_zeros(unique_chains)
    d0dom       = init_chainpairdict_zeros(unique_chains)
    d0dom_max   = init_chainpairdict_zeros(unique_chains)
    d0res       = init_chainpairdict_zeros(unique_chains)
    d0res_max   = init_chainpairdict_zeros(unique_chains)
    d0res_byres = init_chainpairdict_npzeros(unique_chains, numres)

    valid_pair_counts           = init_chainpairdict_zeros(unique_chains)
    dist_valid_pair_counts      = init_chainpairdict_zeros(unique_chains)
    unique_residues_chain1      = init_chainpairdict_set(unique_chains)
    unique_residues_chain2      = init_chainpairdict_set(unique_chains)
    dist_unique_residues_chain1 = init_chainpairdict_set(unique_chains)
    dist_unique_residues_chain2 = init_chainpairdict_set(unique_chains)

    # calculate ipTM/ipSAE with and without PAE cutoff

    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1 == chain2:
                continue

            n0chn[chain1][chain2]=np.sum( chains==chain1) + np.sum(chains==chain2) # total number of residues in chain1 and chain2
            d0chn[chain1][chain2]=calc_d0(n0chn[chain1][chain2], chain_pair_type[chain1][chain2])
            ptm_matrix_d0chn=np.zeros((numres,numres))
            ptm_matrix_d0chn=ptm_func_vec(pae_matrix,d0chn[chain1][chain2])

            valid_pairs_iptm = (chains == chain2)
            valid_pairs_matrix = np.outer(chains == chain1, chains == chain2) & (pae_matrix < pae_cutoff)

            for i in range(numres):


                if chains[i] != chain1:
                    continue

                valid_pairs_ipsae = valid_pairs_matrix[i]  # row for residue i of chain1
                iptm_d0chn_byres[chain1][chain2][i] =  ptm_matrix_d0chn[i, valid_pairs_iptm].mean() if valid_pairs_iptm.any() else 0.0
                ipsae_d0chn_byres[chain1][chain2][i] = ptm_matrix_d0chn[i, valid_pairs_ipsae].mean() if valid_pairs_ipsae.any() else 0.0

                # Track unique residues contributing to the IPSAE for chain1,chain2
                valid_pair_counts[chain1][chain2] += np.sum(valid_pairs_ipsae)
                if valid_pairs_ipsae.any():
                    iresnum=residues[i]['resnum']
                    unique_residues_chain1[chain1][chain2].add(iresnum)
                    for j in np.where(valid_pairs_ipsae)[0]:
                        jresnum=residues[j]['resnum']
                        unique_residues_chain2[chain1][chain2].add(jresnum)

                # Track unique residues contributing to iptm in interface
                valid_pairs = (chains == chain2) & (pae_matrix[i] < pae_cutoff) & (distances[i] < dist_cutoff)
                dist_valid_pair_counts[chain1][chain2] += np.sum(valid_pairs)

                # Track unique residues contributing to the IPTM
                if valid_pairs.any():
                    iresnum=residues[i]['resnum']
                    dist_unique_residues_chain1[chain1][chain2].add(iresnum)
                    for j in np.where(valid_pairs)[0]:
                        jresnum=residues[j]['resnum']
                        dist_unique_residues_chain2[chain1][chain2].add(jresnum)

    #OUT2.write("i   AlignChn ScoredChain  AlignResNum  AlignResType  AlignRespLDDT      n0chn  n0dom  n0res    d0chn     d0dom     d0res   ipTM_pae  ipSAE_d0chn ipSAE_d0dom    ipSAE \n")
    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1 == chain2:
                continue
            residues_1 = len(unique_residues_chain1[chain1][chain2])
            residues_2 = len(unique_residues_chain2[chain1][chain2])
            n0dom[chain1][chain2] = residues_1+residues_2
            d0dom[chain1][chain2] = calc_d0(n0dom[chain1][chain2], chain_pair_type[chain1][chain2])

            ptm_matrix_d0dom = np.zeros((numres,numres))
            ptm_matrix_d0dom = ptm_func_vec(pae_matrix,d0dom[chain1][chain2])

            valid_pairs_matrix = np.outer(chains == chain1, chains == chain2) & (pae_matrix < pae_cutoff)

            # Assuming valid_pairs_matrix is already defined
            n0res_byres_all = np.sum(valid_pairs_matrix, axis=1)
            d0res_byres_all = calc_d0_array(n0res_byres_all, chain_pair_type[chain1][chain2])

            n0res_byres[chain1][chain2] = n0res_byres_all
            d0res_byres[chain1][chain2] = d0res_byres_all

            for i in range(numres):
                if chains[i] != chain1:
                    continue
                valid_pairs = valid_pairs_matrix[i]
                ipsae_d0dom_byres[chain1][chain2][i] = ptm_matrix_d0dom[i, valid_pairs].mean() if valid_pairs.any() else 0.0

                ptm_row_d0res=np.zeros((numres))
                ptm_row_d0res=ptm_func_vec(pae_matrix[i], d0res_byres[chain1][chain2][i])
                ipsae_d0res_byres[chain1][chain2][i] = ptm_row_d0res[valid_pairs].mean() if valid_pairs.any() else 0.0

                # outstring = f'{i+1:<4d}    ' + (
                #     f'{chain1:4}      '
                #     f'{chain2:4}      '
                #     f'{residues[i]["resnum"]:4d}           '
                #     f'{residues[i]["res"]:3}        '
                #     f'{plddt[i]:8.2f}         '
                #     f'{int(n0chn[chain1][chain2]):5d}  '
                #     f'{int(n0dom[chain1][chain2]):5d}  '
                #     f'{int(n0res_byres[chain1][chain2][i]):5d}  '
                #     f'{d0chn[chain1][chain2]:8.3f}  '
                #     f'{d0dom[chain1][chain2]:8.3f}  '
                #     f'{d0res_byres[chain1][chain2][i]:8.3f}   '
                #     f'{iptm_d0chn_byres[chain1][chain2][i]:8.4f}    '
                #     f'{ipsae_d0chn_byres[chain1][chain2][i]:8.4f}    '
                #     f'{ipsae_d0dom_byres[chain1][chain2][i]:8.4f}    '
                #     f'{ipsae_d0res_byres[chain1][chain2][i]:8.4f}\n'
                # )
                #OUT2.write(outstring)

    # Compute interchain ipTM and ipSAE for each chain pair
    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1 == chain2:
                continue

            interchain_values = iptm_d0chn_byres[chain1][chain2]
            max_index = np.argmax(interchain_values)
            iptm_d0chn_asym[chain1][chain2] = interchain_values[max_index]
            iptm_d0chn_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"

            interchain_values = ipsae_d0chn_byres[chain1][chain2]
            max_index = np.argmax(interchain_values)
            ipsae_d0chn_asym[chain1][chain2] = interchain_values[max_index]
            ipsae_d0chn_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"

            interchain_values = ipsae_d0dom_byres[chain1][chain2]
            max_index = np.argmax(interchain_values)
            ipsae_d0dom_asym[chain1][chain2] = interchain_values[max_index]
            ipsae_d0dom_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"

            interchain_values = ipsae_d0res_byres[chain1][chain2]
            max_index = np.argmax(interchain_values)
            ipsae_d0res_asym[chain1][chain2] = interchain_values[max_index]
            ipsae_d0res_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"
            n0res[chain1][chain2]=n0res_byres[chain1][chain2][max_index]
            d0res[chain1][chain2]=d0res_byres[chain1][chain2][max_index]

            # pick maximum value for each chain pair for each iptm/ipsae type
            if chain1 > chain2:
                maxvalue=max(iptm_d0chn_asym[chain1][chain2], iptm_d0chn_asym[chain2][chain1])
                if maxvalue==iptm_d0chn_asym[chain1][chain2]: maxres=iptm_d0chn_asymres[chain1][chain2]
                else: maxres=iptm_d0chn_asymres[chain2][chain1]
                iptm_d0chn_max[chain1][chain2]=maxvalue
                iptm_d0chn_maxres[chain1][chain2]=maxres
                iptm_d0chn_max[chain2][chain1]=maxvalue
                iptm_d0chn_maxres[chain2][chain1]=maxres

                maxvalue=max(ipsae_d0chn_asym[chain1][chain2], ipsae_d0chn_asym[chain2][chain1])
                if maxvalue==ipsae_d0chn_asym[chain1][chain2]: maxres=ipsae_d0chn_asymres[chain1][chain2]
                else: maxres=ipsae_d0chn_asymres[chain2][chain1]
                ipsae_d0chn_max[chain1][chain2]=maxvalue
                ipsae_d0chn_maxres[chain1][chain2]=maxres
                ipsae_d0chn_max[chain2][chain1]=maxvalue
                ipsae_d0chn_maxres[chain2][chain1]=maxres

                maxvalue=max(ipsae_d0dom_asym[chain1][chain2], ipsae_d0dom_asym[chain2][chain1])
                if maxvalue==ipsae_d0dom_asym[chain1][chain2]:
                    maxres=ipsae_d0dom_asymres[chain1][chain2]
                    maxn0=n0dom[chain1][chain2]
                    maxd0=d0dom[chain1][chain2]
                else:
                    maxres=ipsae_d0dom_asymres[chain2][chain1]
                    maxn0=n0dom[chain2][chain1]
                    maxd0=d0dom[chain2][chain1]
                ipsae_d0dom_max[chain1][chain2]=maxvalue
                ipsae_d0dom_maxres[chain1][chain2]=maxres
                ipsae_d0dom_max[chain2][chain1]=maxvalue
                ipsae_d0dom_maxres[chain2][chain1]=maxres
                n0dom_max[chain1][chain2]=maxn0
                n0dom_max[chain2][chain1]=maxn0
                d0dom_max[chain1][chain2]=maxd0
                d0dom_max[chain2][chain1]=maxd0

                maxvalue=max(ipsae_d0res_asym[chain1][chain2], ipsae_d0res_asym[chain2][chain1])
                if maxvalue==ipsae_d0res_asym[chain1][chain2]:
                    maxres=ipsae_d0res_asymres[chain1][chain2]
                    maxn0=n0res[chain1][chain2]
                    maxd0=d0res[chain1][chain2]
                else:
                    maxres=ipsae_d0res_asymres[chain2][chain1]
                    maxn0=n0res[chain2][chain1]
                    maxd0=d0res[chain2][chain1]
                ipsae_d0res_max[chain1][chain2]=maxvalue
                ipsae_d0res_maxres[chain1][chain2]=maxres
                ipsae_d0res_max[chain2][chain1]=maxvalue
                ipsae_d0res_maxres[chain2][chain1]=maxres
                n0res_max[chain1][chain2]=maxn0
                n0res_max[chain2][chain1]=maxn0
                d0res_max[chain1][chain2]=maxd0
                d0res_max[chain2][chain1]=maxd0

    chainpairs=set()
    for chain1 in unique_chains:
        for chain2 in unique_chains:
            if chain1 >= chain2: continue
            chainpairs.add(chain1 + "-" + chain2)

    if write_file:
        path_stem =  f'{str(pdb_path).replace(".cif","")}_{pae_string}_{dist_string}'
        file_path = path_stem + "_new.txt"
        OUT = open(file_path,'w')
        OUT.write("\nChn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model\n")

    rows = [] 

    #PML.write("# Chn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model\n")
    for pair in sorted(chainpairs):
        (chain_a, chain_b) = pair.split("-")
        pair1 = (chain_a, chain_b)
        pair2 = (chain_b, chain_a)
        for pair in (pair1, pair2):
            chain1=pair[0]
            chain2=pair[1]

            residues_1 = len(unique_residues_chain1[chain1][chain2])
            residues_2 = len(unique_residues_chain2[chain1][chain2])
            dist_residues_1 = len(dist_unique_residues_chain1[chain1][chain2])
            dist_residues_2 = len(dist_unique_residues_chain2[chain1][chain2])
            #pairs = valid_pair_counts[chain1][chain2]
            #dist_pairs = dist_valid_pair_counts[chain1][chain2]
            iptm_af = iptm_af3[chain1][chain2]  # symmetric value for each chain pair

            if write_file:
                outstring=f'{chain1}    {chain2}     {pae_string:3}  {dist_string:3}  {"asym":5} ' + (
                    f'{ipsae_d0res_asym[chain1][chain2]:8.6f}    '
                    f'{ipsae_d0chn_asym[chain1][chain2]:8.6f}    '
                    f'{ipsae_d0dom_asym[chain1][chain2]:8.6f}    '
                    f'{iptm_af:5.3f}    '
                    f'{iptm_d0chn_asym[chain1][chain2]:8.6f}    '
                    f'{int(n0res[chain1][chain2]):5d}  '
                    f'{int(n0chn[chain1][chain2]):5d}  '
                    f'{int(n0dom[chain1][chain2]):5d}  '
                    f'{d0res[chain1][chain2]:6.2f}  '
                    f'{d0chn[chain1][chain2]:6.2f}  '
                    f'{d0dom[chain1][chain2]:6.2f}  '
                    f'{residues_1:5d}   '
                    f'{residues_2:5d}   '
                    f'{dist_residues_1:5d}   '
                    f'{dist_residues_2:5d}   \n')
                OUT.write(outstring)
            row = {
                "Chn1": chain1,
                "Chn2": chain2,
                "PAE": int(pae_string),
                "Dist": int(dist_string),
                "Type": "asym",
                "ipSAE": ipsae_d0res_asym[chain1][chain2],
                "ipSAE_d0chn": ipsae_d0chn_asym[chain1][chain2],
                "ipSAE_d0dom": ipsae_d0dom_asym[chain1][chain2],
                "ipTM_af": iptm_af,
                "ipTM_d0chn": iptm_d0chn_asym[chain1][chain2],
                "n0res": int(n0res[chain1][chain2]),
                "n0chn": int(n0chn[chain1][chain2]),
                "n0dom": int(n0dom[chain1][chain2]),
                "d0res": d0res[chain1][chain2],
                "d0chn": d0chn[chain1][chain2],
                "d0dom": d0dom[chain1][chain2],
                "nres1": residues_1,
                "nres2": residues_2,
                "dist1": dist_residues_1,
                "dist2": dist_residues_2,
            }
            rows.append(row)
            #PML.write("# " + outstring)
            if chain1 > chain2:
                residues_1 = max(len(unique_residues_chain2[chain1][chain2]), len(unique_residues_chain1[chain2][chain1]))
                residues_2 = max(len(unique_residues_chain1[chain1][chain2]), len(unique_residues_chain2[chain2][chain1]))
                dist_residues_1 = max(len(dist_unique_residues_chain2[chain1][chain2]), len(dist_unique_residues_chain1[chain2][chain1]))
                dist_residues_2 = max(len(dist_unique_residues_chain1[chain1][chain2]), len(dist_unique_residues_chain2[chain2][chain1]))

                iptm_af_value=iptm_af

                if write_file:
                    outstring=f'{chain2}    {chain1}     {pae_string:3}  {dist_string:3}  {"max":5} ' + (
                        f'{ipsae_d0res_max[chain1][chain2]:8.6f}    '
                        f'{ipsae_d0chn_max[chain1][chain2]:8.6f}    '
                        f'{ipsae_d0dom_max[chain1][chain2]:8.6f}    '
                        f'{iptm_af_value:5.3f}    '
                        f'{iptm_d0chn_max[chain1][chain2]:8.6f}    '
                        f'{int(n0res_max[chain1][chain2]):5d}  '
                        f'{int(n0chn[chain1][chain2]):5d}  '
                        f'{int(n0dom_max[chain1][chain2]):5d}  '
                        f'{d0res_max[chain1][chain2]:6.2f}  '
                        f'{d0chn[chain1][chain2]:6.2f}  '
                        f'{d0dom_max[chain1][chain2]:6.2f}  '
                        f'{residues_1:5d}   '
                        f'{residues_2:5d}   '
                        f'{dist_residues_1:5d}   '
                        f'{dist_residues_2:5d}   \n')
                    OUT.write(outstring)

                row = {
                    "Chn1": chain1,
                    "Chn2": chain2,
                    "PAE": int(pae_string),
                    "Dist": int(dist_string),
                    "Type": "asym",
                    "ipSAE": ipsae_d0res_asym[chain1][chain2],
                    "ipSAE_d0chn": ipsae_d0chn_asym[chain1][chain2],
                    "ipSAE_d0dom": ipsae_d0dom_asym[chain1][chain2],
                    "ipTM_af": iptm_af,
                    "ipTM_d0chn": iptm_d0chn_asym[chain1][chain2],
                    "n0res": int(n0res[chain1][chain2]),
                    "n0chn": int(n0chn[chain1][chain2]),
                    "n0dom": int(n0dom[chain1][chain2]),
                    "d0res": d0res[chain1][chain2],
                    "d0chn": d0chn[chain1][chain2],
                    "d0dom": d0dom[chain1][chain2],
                    "nres1": residues_1,
                    "nres2": residues_2,
                    "dist1": dist_residues_1,
                    "dist2": dist_residues_2,
                }
                rows.append(row)

                #PML.write("# " + outstring)

            #chain_pair= f'color_{chain1}_{chain2}'
            #chain1_residues = f'chain  {chain1} and resi {contiguous_ranges(unique_residues_chain1[chain1][chain2])}'
            #chain2_residues = f'chain  {chain2} and resi {contiguous_ranges(unique_residues_chain2[chain1][chain2])}'
            #PML.write(f'alias {chain_pair}, color gray80, all; color {color1}, {chain1_residues}; color {color2}, {chain2_residues}\n\n')
        
        #OUT.write("\n")
    if write_file:
        OUT.write("\n")
        OUT.close()
        return pd.DataFrame(rows), file_path

    return pd.DataFrame(rows)
