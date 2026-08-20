#!/bin/bash
#SBATCH --job-name=af3_metrics
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=150:00:00
#Fallback path for standalone submission; workflow.sh overrides these via sbatch --output/--error
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

#Strict execution mode
set -Eeuo pipefail

#Usage: sbatch --array=0-<num_splits-1>%<concurrency> collect_metrics_slurm_parallel.sh <folder_path> <suffix> <num_splits> <logs_path>
#CONDA_SH must be set in the environment (see configs/env.cfg.example); sbatch inherits it
#from the submitting shell, which sources configs/env.cfg before calling sbatch.
: "${CONDA_SH:?CONDA_SH is not set — define it in configs/env.cfg (see configs/env.cfg.example)}"
source "$CONDA_SH"
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
