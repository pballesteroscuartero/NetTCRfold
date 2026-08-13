#!/bin/bash
#SBATCH --job-name=af3_metrics
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=150:00:00
#SBATCH --output=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/slurm_%A_%a.out
#SBATCH --error=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/slurm_%A_%a.err

#Strict execution mode
set -Eeuo pipefail

#Usage: sbatch --array=0-<num_splits-1>%<concurrency> collect_metrics_slurm_parallel.sh <folder_path> <suffix> <num_splits> <logs_path>
source /home/projects2/pbacu/utils/Miniconda/etc/profile.d/conda.sh
conda activate structureTCR

folder_path=$1
suffix=$2
num_splits=$3
logs_path=$4

mkdir -p "$logs_path"
exec >"${logs_path}/split${SLURM_ARRAY_TASK_ID}.out" 2>"${logs_path}/split${SLURM_ARRAY_TASK_ID}.err"

python -m structureTCR.metrics.collect_af3metrics_extended_parallel \
    -i "$folder_path" \
    -s "$suffix" \
    --split_idx "$SLURM_ARRAY_TASK_ID" \
    --num_splits "$num_splits"
