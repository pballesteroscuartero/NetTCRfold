#!/bin/bash
#SBATCH --job-name=runAF3_fullworkflow
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500MB
#SBATCH --time=96:00:00
#SBATCH --output=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/slurm_%A_%a.out
#SBATCH --error=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/slurm_%A_%a.err

#Strict execution mode
set -Eeuo pipefail

#Usage example sbatch dataGenerationWorklow.sh + Write any change needed in config file
source /home/projects2/pbacu/utils/Miniconda/etc/profile.d/conda.sh
conda activate structurePipeline

source configs/af3Benchmark_allfolders_5seeds_50diffusion.cfg

#Define paths
suffix_output_inference=$SUFFIX_OUTPUT

output_base=$OUTPUT_DIR
output_savedata="${output_base}/data/af3_output"
input_basejson="${output_base}/data/benchmark_jsonFiles"
output_customjson="${input_basejson}/customJSON/"
output_datageneration="${output_savedata}/small_db/dataPipelineOut/"
output_inference="${output_savedata}/small_db/structInference${suffix_output_inference}/"
output_nettcrstruct_datareformatting="${output_base}/data/nettcrstruc${suffix_output_inference}"

logs_path="${output_base}/logs" 
logs_path_datageneration="${logs_path}/af3_datageneration_workflow/smalldb"      
logs_path_inference="${logs_path}/af3_inference/smalldb"
logs_path_nettcrstruct="${logs_path}/nettcrstruc/"

##Redirect log
name_log="${output_base%/}"
name_log="${name_log##*/}"
exec >"/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/${name_log}${suffix_output_inference}${SLURM_JOB_ID}.out" 2>"/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/${name_log}${suffix_output_inference}${SLURM_JOB_ID}.err"

##1.Perform data preprocessing
if $RUN_JSON_WITH_MSA_TEMPLATE_GENERATION; then
    echo "Performing data preprocessing step: Generating JSON for data generation step with MSA and Template"
    cmd=(
    python scripts/dataPreprocessing.py
    -i "$INPUT_DB"
    -o "$output_base/data"
    )

    if [[ -n "${PARTITION:-}" ]]; then
        cmd+=(-p "$PARTITION")
    fi

    "${cmd[@]}"
    echo "Data preprocessing finished"
else
    echo "Skipping data preprocessing step"
fi

if $RUN_DATA_GENERATION_PIPELINE; then
    combinations=(
    "paired standard"
    "unpaired standard"
    "unpaired onquery"
    )  

    for combo in "${combinations[@]}"; do
        read -r uniprot_msa template_selection_method <<< "$combo"
        echo "Running AF data generation step with MSA (uniprot in ${uniprot_msa}) and Template with selection method ${template_selection_method}:"  
        output_folder="${output_datageneration}/uniprotOn_${uniprot_msa}_template_${template_selection_method}"
        max_array=1000 
  
        for start in $(seq $GLOBAL_START $max_array $TOTAL_TASKS); do
            end=$((${start} + ${max_array} - 1))
            [ $end -gt $TOTAL_TASKS ] && end=$TOTAL_TASKS

            array_length=$((${end} - ${start} + 1))
            echo "Submitting array: $start-$end%$CONCURRENT (length $array_length) for combo: $combo"

            jobid=$(sbatch --parsable --wait \
            --array=1-${array_length}%${CONCURRENT} \
            scripts/runAF3_dataGeneration.sh \
            $ARRAY_MAP $input_basejson/json_msa_template $output_folder $logs_path_datageneration \
            $uniprot_msa $template_selection_method $start)

            echo "Job $jobid finished."
        done

    done
else 
    echo "Skipping AF3 data generation step"
fi

if $RUN_CUSTOM_JSON_GENERATION; then
    echo "Generating custom JSON input files with different MSA and Template settings"
    python scripts/create_custom_json.py  -i "${output_datageneration}" -o $output_customjson
    echo "Custom JSON input generation finished. Files saved in $output_customjson"
else
    echo "Skipping custom JSON input generation"
fi

#2. Perform af3 inference step in selected folders
if $RUN_AF3_INFERENCE; then
    echo "Running AF3 inference step"
    
    if [[ -z "${FOLDERS_INFERENCE:-}" ]]; then
        echo "FOLDERS not provided in config — processing all subfolders in $output_customjson"
        folders=$(find "$output_customjson" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')
        echo "Processing folders: $folders"
    else
        folders=$FOLDERS_INFERENCE
        echo "Using folders from config: $folders"
    fi

    for folder in $folders; do
        folder_path="$output_customjson/$folder"
        echo "Submitting job for $folder"
        jobid=$(sbatch --parsable --wait scripts/runAF3inference.sh "$folder_path" "$output_inference" "$logs_path_inference" "$ARRAY_MAP" "$NUM_SEEDS" "$NUM_DIFFUSION")
        echo "Job $jobid for $folder finished."
    done

    echo "All folders processed sequentially."
else 
    echo "Skipping AF3 inference step"
fi

##3. Collect metrics

if $RUN_METRICS_COLLECTION; then

    echo "Collecting metrics step"
    if [[ -z "${FOLDERS_METRIC_COLLECTION:-}" ]]; then
        echo "FOLDERS not provided in config — processing all subfolders in $output_inference"
        folders_metric_collection=$(find "$output_inference" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')
        echo "Folders to process: $folders_metric_collection"
    else
        folders_metric_collection=$FOLDERS_METRIC_COLLECTION
        echo "Using folders from config. Folders to process: $folders_metric_collection"
    fi

    if $COMPUTE_DOCKQ; then
        echo "Computing DockQ scores"
        for folder in $folders_metric_collection; do
            folder_path="$output_inference/$folder"
            echo "Computing Dockq for files in $folder"
            python scripts/computeDockq.py -i "$folder_path" -t "$TEMPLATE_PATH" -n1 D E -m1 D E -n2 C -m2 C
            echo "Computation for $folder finished."
        done
        echo "DockQ collected for all folders"
    fi

    echo "Collecting all metrics for AF3 generated structures"
    for folder in $folders_metric_collection; do
        folder_path="$output_inference/$folder"
        echo "Collecting metrics for files in $folder"
        python scripts/collect_af3metrics_extended.py -i "$folder_path"
        echo "Computation for $folder finished."

    done

    echo "AF3 metrics collection finished for all folders"
else
    echo "Skipping metrics collection"
    #echo "Metrics collection finished. Results saved in $AF3_OUTPUT_FOLDER/small_db/metrics/"
fi


if $COMPUTE_NETTCRSTRUC; then
    echo "Computing reranking with nettcrstruct"

    if [[ -z "${FOLDERS_NETTCRSTRUC:-}" ]]; then
        echo "FOLDERS not provided in config — processing all subfolders in $output_inference"
        #folders_nettcrstruc=$(find "$output_inference" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')
        mapfile -t folders_nettcrstruc < <(
        find "$output_inference" -mindepth 1 -maxdepth 1 -type d -printf '%f\n'
        )
        folders_nettcrstruc_str="${folders_nettcrstruc[*]}"

        echo "Folders to process: $folders_nettcrstruc"
    else
        #folders_nettcrstruc=$FOLDERS_NETTCRSTRUC
        read -r -a folders_nettcrstruc <<< "$FOLDERS_NETTCRSTRUC"
        folders_nettcrstruc_str="${folders_nettcrstruc[*]}"
        echo "Using folders from config. Folders to process: $folders_nettcrstruc"
    fi

    if $PREPARE_NETTCRSTRUCT_INPUT; then
        for folder in $folders_nettcrstruc; do
            folder_path="$output_inference/$folder"
            echo "Submitting nettcrstruct prepare input job for $folder"
            python scripts/prepare_input_nettcrstruc.py -i $folder_path -o $output_nettcrstruct_datareformatting
        done
    else
        echo "Skipping input reformatting for nettcrstruct"
    fi

    num_folders=${#folders_nettcrstruc[@]}
    echo "Submitting netttcrstruct"
    jobid=$(sbatch --parsable --wait --array=1-${num_folders}%4 scripts/run_nettcrstruc.sh $output_nettcrstruct_datareformatting $logs_path_nettcrstruct $ENSEMBLE $suffix_output_inference $folders_nettcrstruc)
    echo "Job $jobid for nettcrstruct finished."
    python scripts/collect_nettcrstruct_reranking.py -i $output_nettcrstruct_datareformatting -o $output_inference

else
    echo "Skipping nettcrstruc computation"

fi

echo "Pipeline completed"