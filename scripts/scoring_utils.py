import numpy as np
import pandas as pd
from itertools import product, repeat
import biotite.structure as struct

from sequence_utils import get_cdr_indices
from structure_utils import get_sequence_from_chain, extract_cdrs_from_structure, get_interface_atoms


def avg_plddt_interface(structure, atom_plddts, chain1, chain2, cutoff=5.0):
    """
    Compute the average pLDDT for atoms at the interface between chain1 and chain2
    Returns:
        avg_plddt (float): mean pLDDT of interface atoms
    """
    atom_pairs = get_interface_atoms(structure, chain1, chain2, cutoff=cutoff)
    if len(atom_pairs) == 0:
        return np.nan

    interface_atom_ids = np.unique(np.array(atom_pairs).flatten())
    mask_atoms = np.isin(structure.atom_id, interface_atom_ids)
    avg_plddt = np.mean(np.array(atom_plddts)[mask_atoms]) * 0.01

    return avg_plddt

def mean_plddt_perinterface(structure, af):
    chains = np.unique(structure.chain_id)

    # Initialize a square DataFrame to hold average interface pLDDT
    df_plddt = pd.DataFrame(np.nan, index=chains, columns=chains)

    # Iterate over all chain pairs
    for chain1, chain2 in product(chains, chains):
        if chain1 == chain2:
            continue  # skip self-interface
        avg_plddt = avg_plddt_interface(structure, af["atom_plddts"], chain1, chain2, cutoff=5.0)
        df_plddt.loc[chain1, chain2] = avg_plddt

    return df_plddt

def metric_pae_interface(structure, af_output, chain1, chain2, cutoff=5.0):
    """
    Compute residue-level PAE statistics (min and mean) for the interface between chain1 and chain2.
    """
    # Get atom pairs at the interface
    atom_pairs = get_interface_atoms(structure, chain1, chain2, cutoff=cutoff)
    if len(atom_pairs) == 0:
        return np.nan, np.nan
    
    interface_atom_ids1 = [i for i, j in atom_pairs]
    interface_atom_ids2 = [j for i, j in atom_pairs]
    residue_mask_chain1 = struct.get_residue_masks(structure, interface_atom_ids1)
    residue_mask_chain2 = struct.get_residue_masks(structure, interface_atom_ids2)
    residue_indices_chain1 = [structure.res_id[mask][0] for mask in residue_mask_chain1]
    residue_indices_chain2 = [structure.res_id[mask][0] for mask in residue_mask_chain2]
    residue_pairs = list(zip(residue_indices_chain1, residue_indices_chain2))
    unique_pairs = list(dict.fromkeys(residue_pairs))
    
    token_chain_ids = np.array(af_output["token_chain_ids"])
    token_res_ids = np.array(af_output["token_res_ids"])
    pae_matrix = np.array(af_output["pae"])
    pae_values = []
   
    for res1, res2 in unique_pairs:
        # find index of residue in each chain
        idx1 = np.where((token_chain_ids == chain1) & (token_res_ids == res1))[0]
        idx2 = np.where((token_chain_ids == chain2) & (token_res_ids == res2))[0]
        if len(idx1) == 0 or len(idx2) == 0:
            continue  # skip if residue not found
        # There should normally be only one index per residue
        pae_values.append(pae_matrix[idx1[0], idx2[0]])

    if len(pae_values) == 0:
        return np.nan, np.nan
    min_pae = np.min(pae_values)
    mean_pae = np.mean(pae_values)

    return min_pae, mean_pae

def mean_min_pae_perinterface(structure, af):
    chains = np.unique(structure.chain_id)
    df_min_pae = pd.DataFrame(np.nan, index=chains, columns=chains)
    df_mean_pae = pd.DataFrame(np.nan, index=chains, columns=chains)

    # Iterate over all chain pairs
    for chain1, chain2 in product(chains, chains):
        if chain1 == chain2:
            continue  # skip self-interface
        min_pae, mean_pae = metric_pae_interface(structure, af, chain1, chain2, cutoff=5.0)
        df_min_pae.loc[chain1, chain2] = min_pae
        df_mean_pae.loc[chain1, chain2] = mean_pae

    return df_min_pae, df_mean_pae

def compute_plddt_sum(
    plddt: np.ndarray,
    indices: list,
    eps: float = 1e-6,
) -> float:
    """Compute summed PLDDT for a CDR.

    Args:
        plddt (np.ndarray): PLDDT scores.
        indices (list): List of lists with CDR indices.
        eps (float): Epsilon value to prevent division by zero.

    Returns:
        float: Summed PLDDT scores.
    """
    indices = np.concatenate(indices)
    return ((plddt[indices].sum() + eps) / len(indices)) * 0.01


def get_plddt_score(
    structure,
    chain_names: list = ["D", "E", "C", "A"],
    plddt_scores: list = None,
) -> list:
    """Gets pLDDT scores for CDR123ab-peptide residues.
    Args:
        structure: struc.AtomArray with pLDDT attribute.
    
        chain_names: List of chain names in the PDB file in the order TCRa, TCRb, peptide, MHC.


    Returns:
        list: List of PLDDT scores.
    """
    cdrs = extract_cdrs_from_structure(structure)
    structure =  structure[structure.element != 'H']

    #Get plddt scores for peptide
    pep_struct = structure[structure.chain_id == chain_names[2]]
    mask_atoms_pep = np.isin(structure.atom_id, pep_struct.atom_id)
    peptide_plddt = np.array(plddt_scores)[mask_atoms_pep]
    
    #Get plddt scores for TCRa and TCRb
    tra_peptide_indices = np.concatenate(get_cdr_indices(get_sequence_from_chain(structure, chain_names[0]), cdrs[0], cdrs[1], cdrs[2])) #residue index within TCRa chain
    trb_peptide_indices = np.concatenate(get_cdr_indices(get_sequence_from_chain(structure, chain_names[1]), cdrs[3], cdrs[4], cdrs[5]))

    tra_struct = structure[(structure.chain_id == chain_names[0]) & np.isin(structure.res_id, tra_peptide_indices)]
    trb_struct = structure[(structure.chain_id == chain_names[1]) & np.isin(structure.res_id, trb_peptide_indices)]
    mask_atoms_tra = np.isin(structure.atom_id, tra_struct.atom_id)
    tra_plddt = np.array(plddt_scores)[mask_atoms_tra]
    mask_atoms_trb = np.isin(structure.atom_id, trb_struct.atom_id)
    trb_plddt = np.array(plddt_scores)[mask_atoms_trb]
    plddt_cdr_pep = np.concatenate([tra_plddt, trb_plddt, peptide_plddt])
    plddt_mean = np.mean(plddt_cdr_pep) * 0.01
    return plddt_mean 
