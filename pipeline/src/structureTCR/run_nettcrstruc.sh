#!/bin/bash
#SBATCH --job-name=nettcrstruc
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=shard:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#Fallback path for standalone submission; workflow.sh overrides these via sbatch --output/--error
#SBATCH --output=logs/nettcrstruc/slurm_%A_%a.out
#SBATCH --error=logs/nettcrstruc/slurm_%A_%a.err

#Usage example: sbatch run_nettcrstruc.sh ../../data/nettcrstruc/rerank_input/inputModels

input_path=$1 
logs_path=$2
ensemble=$3
log_suffix=$4
folders_str=$5

mapfile -t folders <<< "$folders_str"
folder="${folders[$SLURM_ARRAY_TASK_ID-1]}"

mkdir -p $logs_path
exec >"${logs_path}/${SLURM_ARRAY_TASK_ID}_${folder}_${log_suffix}.out" 2> "${logs_path}/${SLURM_ARRAY_TASK_ID}_${folder}_${log_suffix}.err"

input_folder="${input_path}/${folder}" #It should be a directory containing directories. Each of the param_combs
features_save="${input_path}/features/${folder}"
#features_load_gvp_if1="${input_path}/features/${folder}/gvp_if1_embeddings"
#features_load_gvp="${input_path}/features/${folder}/gvp"

if [[ "$ensemble" == "benchmark" ]]; then
     ensemble_gvp_if1=ensemble_benchmark_gvp_if1_ens
     ensemble_gvp=ensemble_benchmark_gvp_ens
     echo "Using benchmark ensemble: $ensemble_gvp_if1 $ensemble_gvp"
else
     ensemble_gvp_if1=ensemble_binding_gvp_if1_ens
     ensemble_gvp=ensemble_binding_gvp_ens
     echo "Using binding ensemble: $ensemble_gvp_if1 $ensemble_gvp"
fi

bash src/structureTCR/run_nettcrstruc_container.sh $input_folder $features_save $ensemble_gvp_if1 $ensemble_gvp

