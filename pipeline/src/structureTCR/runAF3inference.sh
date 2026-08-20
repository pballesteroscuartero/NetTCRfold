#!/bin/bash
#SBATCH --job-name=runAF3inference
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=shard:2
#SBATCH --cpus-per-task=1
#SBATCH --mem=35g
#SBATCH --time=48:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#SBATCH --output=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/af3_inference/smalldb/slurm_%A_%a.out
#SBATCH --error=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/af3_inference/smalldb/slurm_%A_%a.err

export folder_path=$1
export output_dir=$2
export logs_path=$3
export config=$4
export NUM_SEEDS=$5
export NUM_DIFFUSION=$6
export start_id=$7

folder_name="${folder_path##*/}"
mkdir -p "${logs_path}/${folder_name}"

GLOBAL_TASK_ID=$((${start_id} + ${SLURM_ARRAY_TASK_ID} - 1))
sample=$(awk -v ArrayTaskID=$GLOBAL_TASK_ID '$1==ArrayTaskID {print $2}' $config)
exec >"${logs_path}/${folder_name}/${GLOBAL_TASK_ID}_${sample}_${NUM_SEEDS}seeds_${NUM_DIFFUSION}diffusion.out" 2>"${logs_path}/${folder_name}/${GLOBAL_TASK_ID}_${sample}_${NUM_SEEDS}seeds_${NUM_DIFFUSION}diffusion.err"

export input="${folder_path}/${sample}"
export output="${output_dir}/${folder_name}"
mkdir -p "$output"

bash src/structureTCR/run_alphafold3_tcrpmhcdatabase_tcrdiversity.sh  "$input" "$output" FALSE TRUE "standard" "onquery" "1.1" ${NUM_SEEDS} ${NUM_DIFFUSION}




