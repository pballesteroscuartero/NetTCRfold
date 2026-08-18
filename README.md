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

## Pipeline steps:

1.  RUN_JSON_WITH_MSA_TEMPLATE_GENERATION: This step prepares the input csv in a format for input in the pipeline's data generation step. It decomposes the input datapoints per chain, so the same chains are not processed twice.
    Code: dataPreprocessing.py.
    Inputs:
        - i: Path to the CSV file. The CSV should be in the format (see above).
        -o: Path to save the preprocessed files
        -m: MHC database. Will be used to infer allele sequence if sequence not provided under MHCA_aa
    Outputs:
        - {input_file}_hla_withid: Input csv with an extra assigned HLA column and unique identifiers for each unique chain in the database.
        - chainid_to_array.txt: Mapping of unique chain IDs to slurm array indices. This will be used to process data with slurm arrays.
        - samplename_to_array.txt: Mapping of unique pdb_ids to slurm array indices. This will be used to process data with slurm arrays.
        - jsonFiles: Folder where AF3 inputs will be stored.
            - json_msa_template: Folder containing the inputs for the data generation step. One folder per unique chain, containing the information in AF3's required format.
2. 
