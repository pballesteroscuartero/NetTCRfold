#!/bin/bash
#SBATCH --job-name=runAF3_datageneration
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --time=48:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#SBATCH --output=/home/projects2/pbacu/projects/NetTCRfold/logs/af3_datageneration_workflow/smalldb/slurm_%A_%a.out
#SBATCH --error=/home/projects2/pbacu/projects/NetTCRfold/logs/af3_datageneration_workflow/smalldb/slurm_%A_%a.err

#For afbenchmark array 25, for schumacher array 603, for nettcr array 13956
#Usage example: sbatch runAF3_dataGeneration.sh 

config=$1
json_path=$2
output_dir=$3
logs_path=$4
template_selection_method=$5
start_id=$6

mkdir -p ${logs_path}
GLOBAL_TASK_ID=$((${start_id} + ${SLURM_ARRAY_TASK_ID} - 1))
sample=$(awk -v ArrayTaskID=$GLOBAL_TASK_ID '$1==ArrayTaskID {print $2}' $config)
exec >"${logs_path}/${GLOBAL_TASK_ID}_${sample}_${template_selection_method}.out" 2>"${logs_path}/${GLOBAL_TASK_ID}_${sample}_${template_selection_method}.err"

export input="${json_path}/${sample}"
export output="${output_dir}"
mkdir -p "$output"
bash src/structureTCR/run_alphafold3_tcrpmhcdatabase_tcrdiversity.sh  "$input" "$output" TRUE FALSE "$template_selection_method"
 

