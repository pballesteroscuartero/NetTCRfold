## Installation

1. Clone the repository: 

```
git clone git@github.com:pballesteroscuartero/NetTCRstruc2.git
cd NetTCRstruc2
```

2. Request access to AF3 weights from XXXX and place them inside alphafold3_tcrpmhc folder

3. Download the image for running the modified AF3 pipeline from https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/image link

```
wget -P apptainer \
    https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/image/
```

4. Download the chemical components from https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/chemicalComponents and place them in the following path:

```
wget -P alphafold3_tcrpmhc/src/alphafold3/constants/converters \
  https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/chemicalComponents/ccd.pickle \
  https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/chemicalComponents/chemical_component_sets.pickle
```

5. Download the curated databases with: TODO: TEST IT WORKS AS IT SHOULD

```
wget -r -np -nH --cut-dirs=3 -R "index.html*" -P alphafold3_tcrpmhc/ \
  https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/tcrpmhc_databases/

```

6. Install the environment for running the pipeline:

```
conda env create -f pipeline/NetTCRfold.yml
conda activate NetTCRfold
pip install -e .
```

##NOTE: DockQ version 1 is provided. It is cloned from the original repository: https://github.com/wallnerlab/DockQ/tree/v1.0 


##TODO. Add somewhere the errors: Make sure to download apptainer - otherwise the following error will be seen when running data processing or inference in the examples/logs/af3_datageneration_workflow/*.err: FATAL:   While checking container encryption: could not open image /imagePath/alphafold3_tcrpmhc_cuda126-py312.sif: failed to retrieve path for /imagePath/alphafold3_tcrpmhc_cuda126-py312.sif: lstat /imagePath/alphafold3_tcrpmhc_cuda126-py312.sif: no such file or directory 
weights: Make sure to download the weights - otherwise the following error will be seen when running data processing or inference in the examples/logs/af3_datageneration_workflow/*.err: FATAL:   container creation failed: mount hook function failure: mount /pathrepo/alphafold3_tcrpmhc/weights->/mnt/models error: while mounting /pathRepo/weights: mount source /pathRepo/weights doesn't exist
databases: Make sure to download the tcrpmhc_databases - otherwise the following error will be seen when running data processing or inference in the examples/logs/af3_datageneration_workflow/*.err: FATAL:   container creation failed: mount hook function failure: mount /pathRepo/alphafold3_tcrpmhc/tcrpmhc_databases->/pathRepo/alphafold3_tcrpmhc/tcrpmhc_databases: mount source /pathRepo/alphafold3_tcrpmhc/tcrpmhc_databases doesn't exist

## Input data format:

Required a csv with the columns:
    - Epitope_aa or peptide containing the peptide that should be modeled
    - TRA_aa: Sequence for the TCRA chain. In our case we use the version of the chain trimmed to the variable domain.
    - TRB_aa: Sequence for the TCRb chain. In our case we use the version of the chain trimmed to the variable domain.
    - MHCA_aa (optional): Sequence for the MHC alpha chain variable domain. If not present, it will try to be infered from the allele column (see below).
    - allele (optional): Allele of the MHC. Used to infer MHCA_aa sequence if not present. Will be looked up in the mhc database.
    - pdb_id (optional): pdb_id of the datapoint. 
    - A1, A2, A3, B1, B2, B3 (optional): Sequence of the CDRs1-3 from chain A and B. Used only for naming purposes.
    - name (optional): Unique identifier of each datapoint. If not present, first pdb_id will try to be used as unique identified. If not present, unique identifier will either be peptide_A1_A2_A3_B1_B2_B2 or peptide_TRAaa_TRBaa.

## Config file and environment file

Say that we have examples here etc and brief explanation of the fields

## Pipeline steps:

1.  RUN_JSON_WITH_MSA_TEMPLATE_GENERATION: This step prepares the input csv in a format for input in the pipeline's data generation step. It decomposes the input datapoints per chain, so the same chains are not processed twice.
    Section in config file:
        DATA_DIR=Base path where the project is contained. The input data should be contained in this folder. The generated data by the pipeline will also appear in subfolders of this base folder.
        INPUT_FILE=Path to the CSV file

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

2. RUN_DATA_GENERATION_PIPELINE: MSA generation and template selection per chain. In this step, both unpaired and paired MSA are computed. The user can choose which to select in the next step. The input is an AF3 formatted json of each chain and the output is the MSA and template for that datapoint.
    Section in config file:
        TEMPLATE_SELECTION_METHOD (optional): Select either onquery or standard to set the template selection method that will take place. If not set onquery will be performed
        SUFFIX_DATAGEN (optional): Suffix to append to the datafolder. Useful if multiple settings for dataPipeline are used. E.g Different databases used. If not set, no suffix is set.
        GLOBAL_START (optional): ID from the chain_id_array to start from. Useful if

    Code: runAF3_dataGeneration.sh
    Inputs (in this order): 
        - config: Path to the chain_id_to_array generated in the previous step
        - json_path:Path to the input json files: json_msa_template generated by the previous step
        - output_dir: Directory where to save the generated data. By default it will be in DATA_DIR/data/af3_output/dataGenerationPipeline. A suffix to append to the name can be assigned setting SUFFIX_DATAGEN.
        - logs_path: Path to save the data generation outputs. By default it's in ../DATA_DIR/logs/af3_datageneration_pipeline/
        - template_selection_method: Which template selection method to be used. Set by TEMPLATE_SELECTION_METHOD in config
        - start_id: Which chain ID to start processing. Set by GLOBAL_START in config. If not provided set to 1 by deault.

3. RUN_CUSTOM_JSON_GENERATION: In this step, the desired MSA and template configuration is applied to the data, and the different chains are pooled back together into the original complexes. The output of this step is the data ready to introduce in the af3_inference pipeline.
    Section in config file: 
        MSA_TEMPLATE_COMBINATIONS (optional): Provide the combination one needs to generate in the form of msaType_templateType. The options for MSA type are: unpaired, paired, full (use both paired and unpaired MSA) or no (use No MSA). The options for template are onquery (use sequence for template search), standard (use full MSA profile for template search), no (no template). Options are provided in a string separated by a space. If no options are provided, unpaired_onquery is performed.
    Code: create_custom_json.py
    Inputs:
        - i: Input folder for the reconstruction. The folder named dataPipelineOut{suffix}.
        - o: General folder to save the reconstructed json files. In the pipeline it's set to DATA_DIR/jsonFiles/customJSON
        -d Path for the file created in step 1 containing the unique ID per datapoint {input_file}_hla_withid
        -c combinations string to reconstruct

4. RUN_AF3_INFERENCE: Step to run AF3 inference on the reconstructed JSON files. It takes a folder containing json files and returns the models for each of the datapoints inside that folder.
    Section in config file:
        FOLDERS_INFERENCE (optional):  Space separated list containing the name of the folders within jsonFiles/customJSON that one wants to model. For example, to model the datapoints with unpaired MSA and onQuery template folders inference needs to be set to json_unpairedMSA_onqueryTemplate. If nothing is provided, all the folders in jsonFiles/customJSON are processed.
        SUFFIX_OUTPUT (optional) : Suffix for the AF3 inference folder. If not provided, suffix is set to ""
        NUM_SEEDS: Number of seeds to use for inference. If not provided, one seed is used
        NUM_DIFFUSION: Number of diffusion samples to use for inference. If not provided, five samples are produced.
    Code: runAF3_inference.sh
    Inputs (in this order):
        folder_path: Path to the folder containing the json files to process
        output_inference: Path to the folder where output will be saved
        logs_path: Path to the folder where the log files will be saved
        ARRAY_MAP_INFERENCE: Path to the samplename_to_array.txt produced in the first step
        NUM_SEEDS: Number of seeds
        NUM_DIFFUSSION: Number of diffusion samples
        start: Element from samplename_to_array.txt where we start processing

5. COMPUTE_DOCKQ : If there are solved structures available, the pipeline allows the user to compute the DockQ between each model and the solved structure, with respect to the TCR and the peptideMHC complex.
    Section in config file:
        TEMPLATE_PATH (required): Points to the folder where the solved-structures templates are stored. The templates should be truncated in the same way as the inputs. And stored under the name: {pdb_id}.trunc.fit.pdb, with pdb_id matching the name by which the models are saved.
    Code: computeDockq.py
    Inputs:
        -i : Input folder path containing all the pdb_ids to process. It's the output of the inference step i.e json_unpairedMSA_onqueryTemplate
        -t : Template path. Matches the TEMPLATE_PATH in the config file
        -d : Path to the dockQ repo. It needs to be version 1.0
        -n1 : Name of the chains to evaluate (TCR) in the native structure
        -m1 : Name of the chains to evaluate (TCR) in the models
        -n2: Name to the chains to evaluate (pMHC) in the native structure
        -m2: Name to the chains to evaluate (pMHC) in the model
        -s: Suffic to append to the file names
    Ouput: Saves in each pdb_id/sample_seed folder the dockQ output

6. RUN_METRICS_COLLECTION: Collects a set of metrics for model selection and target selection
    Section in config file:
        FOLDERS_METRIC_COLLECTION (optional): Folders to collect the metrics for. If not provided, all the folders present in structInference folder are evaluated. If multiple folders should be passed, provide them as a blank space separated string. I.e: "json_pairedMSA_onqueryTemplate json_unpairedMSA_onqueryTemplate"
        NUM_METRIC_SPLITS (optional): Number of parallel SLURM array tasks to split metrics collection into. This helps to speed up collection if many datapoints are present. By default set to one.
        CONCURRENT_METRICS (optional): Maximum concurrent tasks from NUM_METRIC_SPLITS. Set according to your resources. By default set to one.
    Code: collect_af3metrics_extended_parallel.py called via collect_metrics_slurm_parallel.sh + combine_metrics_onefile.py
    Inputs for collect_metrics_slurm_parallel.sh (in this order):
        - folder_path: Path to folder containing the inference. structureInference...
        - suffix: Suffix to append to the metric files
        NUM_METRIC_SPLITS: Number of splits when computing the metrics
        metrics_logs: Path to the log where each of the logs for the splits will be saved

    Outpus: The output of this step is a collection of metrics saved within each folder i.e json_unpairedMSA_onqueryTemplate called "collected_af3metrics.csv" containing all the metrics for each datapoint

    Inputs for collect_af3metrics_extended_parallel.py:
    -i: Path to the specific folder to evaluate i.e json_unpairedMSA_onqueryTemplate
    -s: Suffic to append to the metric files
    --split_idx: Specific split we are processing
    --num_splits: Total number of splits
    Note: DockQ is picked up automatically if present on disk (dockQ_metrics_*.json files next to each model), independently of the -s suffix. This means COMPUTE_DOCKQ and RUN_METRICS_COLLECTION can be run in either order/separately: as long as DockQ has been computed at some point for a folder, the dockq column will be populated when metrics are collected for it, even if this particular run has COMPUTE_DOCKQ=false. It will only be empty if no DockQ has been computed yet for that folder.
    Inputs for combine_metrics_onefile.py: Use this code if multiple folders were processed and output should be collapsed into a single file.
        -i: General folder where inference samples are stored i.e: structureInference
        -s: Suffix to append to the file
    Output: A file in the structureInference folder called allresults_merged.csv containing the metrics for all parameter combinations.

    Code: expandMetrics.py. Run after combine_metrics_onefile.py on its output. The per-datapoint metrics are stored as nested per-chain-pair dictionaries (e.g. chain_pair_iptm, chain_pair_pae_min, cdr_metric_mean_chain, ipsae, ipsae_d0chn, ipsae_d0dom); this step flattens them into one column per chain pair (e.g. TRA_TRB_ipsae) so the results can be filtered/plotted without parsing dictionaries.
    Inputs for expandMetrics.py:
        -i/--input_csv: Path to the combined metrics file to expand. In the pipeline this is allresults_merged{suffix}.csv, produced by combine_metrics_onefile.py
        -o/--output_csv: Path to save the expanded csv
    Output: A file in the structureInference folder called allresults_merged_expanded.csv, with the chain-pair metrics expanded into individual columns instead of nested dictionaries.


