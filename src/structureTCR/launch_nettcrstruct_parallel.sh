#!/bin/bash
#SBATCH --job-name=nettcrstruc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=36:00:00
#SBATCH --output=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/nettcrstruc/slurm_%A_%a.out
#SBATCH --error=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/nettcrstruc/slurm_%A_%a.err

input_path="${1%/}"
logs_path=$2
ensemble=$3
log_suffix=$4
folders_str=$5

num_jobs=10

if [[ -z "$input_path" || -z "$logs_path" || -z "$ensemble" || -z "$log_suffix" || -z "$folders_str" ]]; then
    echo "Usage: $0 <input_path> <logs_path> <ensemble> <log_suffix> <folders_str> [num_jobs]"
    exit 1
fi

mkdir -p "$logs_path"
mapfile -t folders <<< "$folders_str"

echo "Found ${#folders[@]} parent folder(s) to process"
echo ""

all_submitted_jobs=()

for folder in "${folders[@]}"; do
    folder="${folder%"${folder##*[! $'\t'$'\r']}"}"

    if [[ -z "$folder" ]]; then
        continue
    fi

    parent_path="${input_path}/${folder}"

    if [[ ! -d "$parent_path" ]]; then
        echo "WARNING: '$parent_path' does not exist, skipping."
        continue
    fi

    # Collect subfolders within this parent
    mapfile -t all_subfolders < <(find "$parent_path" -mindepth 1 -maxdepth 1 -type d | sort)
    total=${#all_subfolders[@]}

    if [[ $total -eq 0 ]]; then
        echo "WARNING: no subfolders found in '$parent_path', skipping."
        continue
    fi

    echo "Parent: $folder ($total subfolders, splitting across $num_jobs jobs)"

    parent_jobs=()
    for ((job_id=0; job_id<num_jobs; job_id++)); do
    #for ((job_id=2; job_id<=2; job_id++)); do
        batch_subfolders=()
        for ((i=job_id; i<total; i+=num_jobs)); do
            batch_subfolders+=("${all_subfolders[$i]}")
        done

        if [[ ${#batch_subfolders[@]} -eq 0 ]]; then
            continue
        fi

        # Permanent split directory — stays here forever
        split_dir="${parent_path}/job_split_${job_id}"
        mkdir -p "$split_dir"

        for subfolder in "${batch_subfolders[@]}"; do
            mv "$subfolder" "${split_dir}/$(basename "$subfolder")"
        done

        echo "  Job $((job_id+1))/${num_jobs}: ${#batch_subfolders[@]} subfolders → $split_dir"

        #job_id_slurm=$(sbatch --parsable \
        #    scripts/run_nettcrstruc_batch.sh \
        #    "$split_dir" \
        #    "$logs_path" \
        #    "$ensemble" \
        #    "${log_suffix}_${folder}_job${job_id}")

        job_id_slurm=$(sbatch --parsable \
            src/structureTCR/run_nettcrstruc_batch.sh \
            "$split_dir" \
            "$input_path" \
            "$folder" \
            "$job_id" \
            "$logs_path" \
            "$ensemble" \
            "${log_suffix}_${folder}_job${job_id}")
            
        parent_jobs+=("$job_id_slurm")
        all_submitted_jobs+=("$job_id_slurm")
    done

    echo "  Submitted job IDs: ${parent_jobs[*]}"
    echo ""
done

echo "All jobs submitted. Total: ${#all_submitted_jobs[@]} jobs"
echo "Job IDs: ${all_submitted_jobs[*]}"
echo ""
echo "Monitor with:  squeue -u \$USER"