#!/bin/bash
#SBATCH --job-name=runAF3_datageneration
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --time=48:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#Fallback path for standalone submission; workflow.sh overrides these via sbatch --output/--error
#SBATCH --output=logs/af3_datageneration_workflow/slurm_%A_%a.out
#SBATCH --error=logs/af3_datageneration_workflow/slurm_%A_%a.err

config=$1
json_path=$2
output_dir=$3
logs_path=$4
msa_type=$5
template_selection_method=$6
start_id=$7

mkdir -p ${logs_path}
GLOBAL_TASK_ID=$((${start_id} + ${SLURM_ARRAY_TASK_ID} - 1))
sample=$(awk -v ArrayTaskID=$GLOBAL_TASK_ID '$1==ArrayTaskID {print $2}' $config)
exec >"${logs_path}/${GLOBAL_TASK_ID}_${sample}_${msa_type}_${template_selection_method}.out" 2>"${logs_path}/${GLOBAL_TASK_ID}_${sample}_${msa_type}_${template_selection_method}.err"

export input="${json_path}/${sample}"
export output="${output_dir}"
mkdir -p "$output"
bash src/structureTCR/run_alphafold3_tcrpmhcdatabase_tcrdiversity.sh  "$input" "$output" TRUE FALSE "$msa_type" "$template_selection_method"
 

