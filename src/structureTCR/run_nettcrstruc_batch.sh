#!/bin/bash
#SBATCH --job-name=nettcrstruc
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=shard:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=150:00:00
#SBATCH --output=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/nettcrstruc/slurm_%A_%a.out
#SBATCH --error=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/nettcrstruc/slurm_%A_%a.err

# input_path is a temp dir containing symlinks to this job's assigned subfolders.
# The container sees it as a normal parent folder and processes all subfolders in it.

#input_path="${1%/}"
#logs_path=$2
#ensemble=$3
#log_suffix=$4

split_dir="${1%/}"
root_input_path="${2%/}"
parent_folder=$3
job_id=$4
logs_path=$5
ensemble=$6
log_suffix=$7

mkdir -p "$logs_path"

exec > "${logs_path}/${log_suffix}.out" \
     2> "${logs_path}/${log_suffix}.err"

echo "Job ${SLURM_JOB_ID} started at $(date)"
echo "Processing folder: $split_dir"
echo ""

# Resolve ensemble names
if [[ "$ensemble" == "benchmark" ]]; then
    ensemble_gvp_if1=ensemble_benchmark_gvp_if1_ens
    ensemble_gvp=ensemble_benchmark_gvp_ens
    echo "Using benchmark ensemble: $ensemble_gvp_if1 $ensemble_gvp"
else
    ensemble_gvp_if1=ensemble_binding_gvp_if1_ens
    ensemble_gvp=ensemble_binding_gvp_ens
    echo "Using binding ensemble: $ensemble_gvp_if1 $ensemble_gvp"
fi

#features_save="${input_path}/features"
features_save="${root_input_path}/features/${parent_folder}/job_split_${job_id}"
mkdir -p "$features_save"

bash src/structureTCR/run_nettcrstruc_container.sh "$split_dir" "$features_save" "$ensemble_gvp_if1" "$ensemble_gvp"

exit_code=$?
if [[ $exit_code -ne 0 ]]; then
    echo "ERROR: run_nettcrstruc_container.sh failed (exit $exit_code)"
else
    echo "Done: $split_dir"
fi

echo ""
echo "Job ${SLURM_JOB_ID} finished at $(date)"