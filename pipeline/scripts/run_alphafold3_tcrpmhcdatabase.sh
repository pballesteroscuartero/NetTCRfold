#! /bin/bash
# static variables (expect alphafold3 resources to be in directory as this file)
export AF3_RESOURCES_DIR="/home/projects2/pbacu/repositories/alphafold3/"
export AF3_SRC=${AF3_RESOURCES_DIR}
export AF3_IMAGE=${AF3_RESOURCES_DIR}/image/alphafold3_tcrpmhc_cuda126-py312.sif
export AF3_MODEL_PARAMETERS_DIR=${AF3_RESOURCES_DIR}/weights
export AF3_DATABASES_DIR=${AF3_RESOURCES_DIR}/tcrpmhc_databases


# user variables
export AF3_INPUTDIR=$1
export AF3_OUTPUTDIR=$2
export DATA_PIPELINE=$3
export INFERENCE=$4
export msamode="${5:-paired}"
export template_mode="${6:-standard}"
export NUM_SEEDS="${7:-25}"
export NUM_DIFFUSION="${8:-5}"

MSA_ARG=()
TEMPLATE_ARG=()
SEED_ARG=()

if [ "$msamode" = "unpaired" ]; then
    MSA_ARG+=(--unpaired_with_uniprot)
fi

if [ "$template_mode" = "onquery" ]; then
    TEMPLATE_ARG+=(--only_query_for_template)
fi

if [ "$NUM_SEEDS" -gt 1 ]; then
    SEED_ARG+=(--num_seeds="$NUM_SEEDS")
fi

#FULL_DB+=(--small_bfd_database_path=/mnt/tcrpmhc_databases/bfd-first_non_consensus_sequences.fasta --mgnify_database_path=/mnt/tcrpmhc_databases/mgnify.fasta --uniprot_cluster_annot_database_path=/mnt/tcrpmhc_databases/uniprot_all_2021_04.fa --uniref90_database_path=/mnt/tcrpmhc_databases/uniref90_2022_05.fa --seqres_database_path=/mnt/tcrpmhc_databases/pdb_seqres_2022_09_28.fasta)
#MSARED_DB+=(--seqres_database_path=/mnt/tcrpmhc_databases/pdb_seqres_2022_09_28.fasta)
#TEMPLATERED_DB+=(--small_bfd_database_path=/mnt/tcrpmhc_databases/bfd-first_non_consensus_sequences.fasta --mgnify_database_path=/mnt/tcrpmhc_databases/mgnify.fasta --uniprot_cluster_annot_database_path=/mnt/tcrpmhc_databases/uniprot_all_2021_04.fa --uniref90_database_path=/mnt/tcrpmhc_databases/uniref90_2022_05.fa --seqres_database_path=/mnt/tcrpmhc_databases/pdb_seqres_tcrpmhc_iedb_stcrdab.txt)
#ALLRED_DB+=(--seqres_database_path=/mnt/tcrpmhc_databases/pdb_seqres_tcrpmhc_iedb_stcrdab.txt)
#ALLRED_DB_HMM+=(--seqres_database_path=/mnt/tcrpmhc_databases/pdb_seqres_tcrpmhc.txt)
#ALLRED_DB_IEDB_STCRDAB+=(--seqres_database_path=/mnt/tcrpmhc_databases/pdb_seqres_tcrpmhc_iedb_mhc_tcr_stcrdab.txt)

# print paths 
echo AlphaFold3 resource: $AF3_RESOURCES_DIR
echo AlphaFold3 source directory $AF3_SRC
echo AlphaFold3 image: $AF3_IMAGE
echo AlphaFold3 model parameters: $AF3_MODEL_PARAMETERS_DIR
echo AlphaFold3 database: $AF3_DATABASES_DIR
echo AlphaFold3 input: $AF3_INPUTDIR
echo AlphaFold3 output: $AF3_OUTPUTDIR

#run AlphaFold3 apptainer (for debugging use apptainer shell $AF3_IMAGE)
apptainer exec \
     --nv \
     --bind $AF3_INPUTDIR:/mnt/af_input \
     --bind $AF3_OUTPUTDIR:/mnt/af_output \
     --bind $AF3_MODEL_PARAMETERS_DIR:/mnt/models \
     --bind $AF3_DATABASES_DIR:/mnt/tcrpmhc_databases \
     --bind $AF3_SRC:/mnt/af_source \
     $AF3_IMAGE \
     python /mnt/af_source/run_alphafold_tcrpmhc.py \
     --json_path=/mnt/af_input/alphafold_input.json \
     --model_dir=/mnt/models \
     --db_dir=/mnt/tcrpmhc_databases \
     --output_dir=/mnt/af_output \
     --run_data_pipeline=$DATA_PIPELINE \
     --run_inference=$INFERENCE \
     --num_diffusion_samples=$NUM_DIFFUSION \
     "${SEED_ARG[@]}" \
     "${TEMPLATE_ARG[@]}" \
     "${MSA_ARG[@]}"