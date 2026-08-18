## Installation

conda env create -f structureTCR.yml


## Data format:

Required a csv with the columns:
    - Epitope_aa or peptide containing the peptide that should be modeled
    - TRA_aa: Sequence for the TCRA chain. In our case we use the version of the chain trimmed to the variable domain.
    - TRB_aa: Sequence for the TCRb chain. In our case we use the version of the chain trimmed to the variable domain.
    - MHCA_aa (optional): Sequence for the MHC alpha chain variable domain. If not present, it will try to be infered from the allele column (see below).
    - allele (optional): Allele of the MHC. Used to infer MHCA_aa sequence if not present. Will be looked up in the mhc database.
    - pdb_id (optional): pdb_id of the datapoint. 
    - A1, A2, A3, B1, B2, B3 (optional): Sequence of the CDRs1-3 from chain A and B. Used only for naming purposes.
    - name (optional): Unique identifier of each datapoint. If not present, first pdb_id will try to be used as unique identified. If not present, unique identifier will either be peptide_A1_A2_A3_B1_B2_B2 or peptide_TRAaa_TRBaa.
