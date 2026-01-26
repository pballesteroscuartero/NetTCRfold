#!/bin/bash
#SBATCH --job-name=nettcrstruc
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=shard:12
#SBATCH --cpus-per-task=10
#SBATCH --mem=15M
#SBATCH --time=36:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#SBATCH --output=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/nettcrstruc/slurm_%A_%a.out
#SBATCH --error=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/nettcrstruc/slurm_%A_%a.err

#Usage example: sbatch run_nettcrstruc.sh ../../data/nettcrstruc/rerank_input/inputModels
source /home/projects2/pbacu/utils/Miniconda/etc/profile.d/conda.sh
conda activate nettcrstruc

input_path=$1 
logs_path=$2
ensemble=$3
log_suffix=$4
folders_string=$5

read -r -a folders <<< "$folders_string"
folder=${folders[$((${SLURM_ARRAY_TASK_ID}-1))]}

#config=/home/projects2/pbacu/projects/structureTCR/schumacherDataset/code/5runNettcrstruc/folder_to_samplename.txt
#folder=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

exec >"${logs_path}/${SLURM_ARRAY_TASK_ID}_${folder}_${log_suffix}.out" 2> "${logs_path}/${SLURM_ARRAY_TASK_ID}_${folder}_${log_suffix}.err"

input_folder="${input_path}/${folder}" #It should be a directory containing directories. Each of the param_combs
features_save="${input_path}/features/${folder}"
features_load_gvp_if1="${input_path}/features/${folder}/gvp_if1_embeddings"
features_load_gvp="${input_path}/features/${folder}/gvp"

if [[ "$ensemble" == "benchmark" ]]; then
     ensemble_gvp_if1=ensemble_benchmark_gvp_if1_ens
     ensemble_gvp=ensemble_benchmark_gvp_ens
     echo "Using benchmark ensemble: $ensemble_gvp_if1 $ensemble_gvp"
else
     ensemble_gvp_if1=ensemble_binding_gvp_if1_ens
     ensemble_gvp=ensemble_binding_gvp_ens
     echo "Using binding ensemble: $ensemble_gvp_if1 $ensemble_gvp"
fi

echo "Creating features"
python3 /home/projects2/pbacu/repositories/NetTCR-struc/nettcrstruc/scripts/create_geometric_features.py -i $input_folder  -o $features_save  -n 2 -d cuda --chain_names D E C A 

# echo "Reranking with gvp_if1"

python3 /home/projects2/pbacu/repositories/NetTCR-struc/nettcrstruc/scripts/rerank_docking_poses.py input_dir=$input_folder \
     processed_dir=$features_load_gvp_if1 \
     name=ensemble_gvp_if1 \
     ensemble=$ensemble_gvp_if1 \
     chain_names=[D,E,C,A]

echo "Reranking with gvp"
python3 /home/projects2/pbacu/repositories/NetTCR-struc/nettcrstruc/scripts/rerank_docking_poses.py input_dir=$input_folder \
    processed_dir=$features_load_gvp \
    name=ensemble_gvp \
    ensemble=$ensemble_gvp \
    chain_names=[D,E,C,A] \
    +ensemble.model_0.model.node_in_dim=[30,3] \
    +ensemble.model_1.model.node_in_dim=[30,3] \
    +ensemble.model_2.model.node_in_dim=[30,3] \
    +ensemble.model_3.model.node_in_dim=[30,3] \
    +ensemble.model_4.model.node_in_dim=[30,3] \
    +ensemble.model_5.model.node_in_dim=[30,3] \
    +ensemble.model_6.model.node_in_dim=[30,3] \
    +ensemble.model_7.model.node_in_dim=[30,3] \
    +ensemble.model_8.model.node_in_dim=[30,3] \
    +ensemble.model_9.model.node_in_dim=[30,3]



