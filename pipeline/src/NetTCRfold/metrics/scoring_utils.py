import numpy as np
import pandas as pd
import biotite.structure as struct
from NetTCRfold.metrics.structure_utils import extract_cdrs_from_structure, get_interface_atoms

def is_interface(structure):
    """
    Check whether we have interface per chains
    Returns:
        df_interface (DataFrame): Whether there is interface or not
    """

    #Restrict only to D-C and E-C interfaces
    chain_pairs = [("D", "C"), ("E", "C")]
    chains = ["D", "C", "E"]
    df_interface = pd.DataFrame(np.nan, index=chains, columns=chains)

    for chain1, chain2 in chain_pairs:
        atom_pairs = get_interface_atoms(structure, chain1, chain2, cutoff=5)
        if len(atom_pairs) == 0:
            interface = 0
        else:
            interface = 1
        df_interface.loc[chain1, chain2] = interface

    return df_interface

def compute_chain_pair_metrics(pae, token_chain_ids, token_res_ids, chain1, chain2, cdr_residues=None, cdr_chain=None):
    
    compute_cdr = cdr_residues is not None
    mask1 = token_chain_ids == chain1
    mask2 = token_chain_ids == chain2
    if not np.any(mask1) or not np.any(mask2):
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    pae_submatrix = pae[np.ix_(mask1, mask2)]
    cdr_mean = np.nan 
    if compute_cdr:
        if cdr_chain == chain1:
            if np.any(cdr_residues):
                pae_cdr = pae_submatrix[cdr_residues, :]
                cdr_mean = float(np.mean(pae_cdr))
        elif cdr_chain == chain2:

            if np.any(cdr_residues):
                pae_cdr = pae_submatrix[:, cdr_residues]
                cdr_mean = float(np.mean(pae_cdr))
    return cdr_mean

def pae_metrics_perchain_pair(structure, af_arrays, cdrs=None):
    pae = af_arrays["pae"]
    token_chain_ids = af_arrays["token_chain_ids"]
    token_res_ids = af_arrays["token_res_ids"]

    chains = np.unique(structure.chain_id)
    #df_cdr_metric_min = pd.DataFrame(np.nan, index=chains, columns=chains)
    df_cdr_metric_mean = pd.DataFrame(np.nan, index=chains, columns=chains)

    if cdrs is None:
        cdrs = extract_cdrs_from_structure(structure)
    #cdrs = extract_cdrs_from_structure(structure) #This returns the indices of the CDRs in the chains, but they are positional indices within the chain, not residue IDs   
    cdr3_a = cdrs[2]  # CDR3 of TCRa
    cdr3_b = cdrs[5]  # CDR3 of T
    cdr12_a = np.concatenate(cdrs[:2])  # CDR1 and CDR2 of TCRa
    cdr12_b = np.concatenate(cdrs[3:5])  # C

    CDR_MAP = {
        frozenset({"D", "C"}): ("D", cdr3_a),  
        frozenset({"E", "C"}): ("E", cdr3_b),
        frozenset({"D", "A"}): ("D", cdr12_a),
        frozenset({"E", "A"}): ("E", cdr12_b),
    }

    for chain1 in chains:
        for chain2 in chains:
            if chain1 == chain2:
                continue
            pair = frozenset({chain1, chain2})
            cdr_chain, cdr_residues = CDR_MAP.get(pair, (None, None))
            cdr_mean = compute_chain_pair_metrics(pae, token_chain_ids, token_res_ids, chain1, chain2, cdr_residues = cdr_residues, cdr_chain=cdr_chain)
            df_cdr_metric_mean.loc[chain1, chain2] = cdr_mean
    return  df_cdr_metric_mean


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
    plddt_scores: list = None
) -> float:
    """Gets pLDDT scores for CDR123ab-peptide residues.
    Args:
        structure: struc.AtomArray with pLDDT attribute.
    
        chain_names: List of chain names in the PDB file in the order TCRa, TCRb, peptide, MHC.


    Returns:
        float: Mean PLDDT score.
    """
    cdrs = extract_cdrs_from_structure(structure) #modified so that this is the idx

    #Get plddt scores for peptide
    pep_struct = structure[structure.chain_id == chain_names[2]]

    tra_cdrs_indices = np.concatenate(cdrs[:3])  # CDR1, CDR2, CDR3 of TCRa. These are a positional mask within the chain (at the residue level)
    trb_cdrs_indices = np.concatenate(cdrs[3:])  # CDR1, CDR2, CDR3 of TCRb. These are a positional mask within the chain (at the residue level)

    residues_chainE = structure[structure.chain_id == chain_names[1]].res_id
    residues_chainD = structure[structure.chain_id == chain_names[0]].res_id

    residues_keep_chainE = residues_chainE[trb_cdrs_indices]
    residues_keep_chainD = residues_chainD[tra_cdrs_indices]

    chainE =structure[structure.chain_id == chain_names[1]]
    chainD =structure[structure.chain_id == chain_names[0]]

    mask_cdr_tra = np.isin(chainD.res_id, residues_keep_chainD)
    mask_cdr_trb = np.isin(chainE.res_id, residues_keep_chainE)

    cdr_tra_struct = chainD[mask_cdr_tra]  
    cdr_trb_struct = chainE[mask_cdr_trb] 

    cdr_struct = struct.concatenate((cdr_tra_struct, cdr_trb_struct))
    mask_cdr = np.isin(structure.atom_id, cdr_struct.atom_id)
    mask_atoms_pep = np.isin(structure.atom_id, pep_struct.atom_id)
    plddt_cdr_pep = np.concatenate([plddt_scores[mask_cdr], plddt_scores[mask_atoms_pep]])
    plddt_mean_cdr_pep = np.mean(plddt_cdr_pep) * 0.01
    return plddt_mean_cdr_pep


def parse_ipsae(df, metrics=None, type_filter="asym"):
    if metrics is None:
        metrics = ["ipSAE", "ipSAE_d0chn", "ipSAE_d0dom"]
        
    df_asym = df[df["Type"] == type_filter].copy()
    all_chains = sorted(set(df_asym["Chn1"]).union(df_asym["Chn2"]))
    metric_dict = {
        metric: {chain1: {chain2: np.nan for chain2 in all_chains} for chain1 in all_chains}
        for metric in metrics
    }

    for _, row in df_asym.iterrows():
        chain1, chain2 = row["Chn1"], row["Chn2"]
        for metric in metrics:
            metric_dict[metric][chain1][chain2] = row[metric]

    return metric_dict



