import numpy as np
import biotite.structure as struc

from sequence_utils import get_cdr_from_sequence, convert_3to1

def get_interface_atoms(structure, chain1, chain2, cutoff=5.0):
    atoms1 = structure[(structure.chain_id == chain1) & (structure.element != 'H')]
    atoms2 = structure[(structure.chain_id == chain2) & (structure.element != 'H')]
    coords1 = atoms1.coord
    coords2 = atoms2.coord
    dists = struc.distance(coords1[:, np.newaxis], coords2)
    idx1, idx2 = np.where(dists <= cutoff)
    atom_pairs_ids = list(zip(structure.atom_id[idx1], structure.atom_id[idx2]))
    return atom_pairs_ids

def get_sequence_from_chain(structure, chain_id: str) -> str:
    """Get the sequence of a chain in a structure.

    Args:
        structure (AtomArray): The structure.
        chain_id (str): The chain ID.

    Returns:
        str: The sequence of the chain in 1-letter code.
    """
    return "".join(
        struc.apply_residue_wise(
            structure[structure.chain_id == chain_id],
            structure[structure.chain_id == chain_id].res_name,
            convert_3to1,
            axis=None,
        )
    )

def extract_cdrs_from_structure(structure):
    tra_seq = get_sequence_from_chain(structure, "D")
    trb_seq = get_sequence_from_chain(structure, "E")

    cdrs = [
        get_cdr_from_sequence(tra_seq, 1)[0],
        get_cdr_from_sequence(tra_seq, 2)[0],
        get_cdr_from_sequence(tra_seq, 3)[0],
        get_cdr_from_sequence(trb_seq, 1)[0],
        get_cdr_from_sequence(trb_seq, 2)[0],
        get_cdr_from_sequence(trb_seq, 3)[0],
    ]

    return cdrs