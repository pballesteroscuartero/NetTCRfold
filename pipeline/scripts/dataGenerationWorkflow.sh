#!/bin/bash
#SBATCH --job-name=runAF3_fullworkflow
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500MB
#SBATCH --time=96:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#SBATCH --output=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/slurm_%A_%a.out
#SBATCH --error=/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/slurm_%A_%a.err

#Strict execution mode
set -Eeuo pipefail

#Usage example sbatch dataGenerationWorklow.sh + Write any change needed in config file
#source /home/projects2/pbacu/utils/Miniconda/etc/profile.d/conda.sh
#conda activate structurePipeline

source configs/swapped_negatives.cfg

src=/home/projects2/pbacu/projects/structureTCR/structurePipeline
mhcdb=$src/databases/mhc_sequences
#Define paths
suffix_output_inference=$SUFFIX_OUTPUT

output_base=$OUTPUT_DIR
output_savedata="${output_base}/data/af3_output"
input_basejson="${output_base}/data/jsonFiles"
output_customjson="${input_basejson}/customJSON/"
output_datageneration="${output_savedata}/small_db/dataPipelineOut/"
output_inference="${output_savedata}/small_db/structInference${suffix_output_inference}/"
output_nettcrstruct_datareformatting="${output_base}/data/nettcrstruc${suffix_output_inference}"

logs_path="${output_base}/logs" 
logs_path_datageneration="${logs_path}/af3_datageneration_workflow/smalldb"      
logs_path_inference="${logs_path}/af3_inference/smalldb"
logs_path_nettcrstruct="${logs_path}/nettcrstruc"

mkdir -p $output_customjson

##Redirect log
name_log="${output_base%/}"
name_log="${name_log##*/}"
exec >"/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/${name_log}${suffix_output_inference}_${SLURM_JOB_ID}.out" 2>"/home/projects2/pbacu/projects/structureTCR/structurePipeline/logs/${name_log}${suffix_output_inference}_${SLURM_JOB_ID}.err"

##1.Perform data preprocessing
if $RUN_JSON_WITH_MSA_TEMPLATE_GENERATION; then
    echo "Performing data preprocessing step: Generating JSON for data generation step with MSA and Template"

    cmd=(
    python /mnt/source/scripts/dataPreprocessing_optimized.py
    -i /mnt/input_db/$INPUT_FILE
    -o /mnt/output/data
    )
    if [[ -n "${PARTITION:-}" ]]; then
        cmd+=(-p "$PARTITION")
    fi

    apptainer exec \
        --bind "$INPUT_DB:/mnt/input_db" \
        --bind "$output_base:/mnt/output" \
        --bind "$src:/mnt/source" \
        --bind "$mhcdb:/mnt/mhc_sequences" \
        "$IMAGE" \
        "${cmd[@]}"

    echo "Data preprocessing finished"
else
    echo "Skipping data preprocessing step"
fi

if $RUN_DATA_GENERATION_PIPELINE; then
    combinations=(
    "paired onquery"
    "unpaired onquery"
    #"unpaired standard"
    #"nomsa onquery"
    )  
    for combo in "${combinations[@]}"; do
        read -r uniprot_msa template_selection_method <<< "$combo"
        if [[ "$uniprot_msa" == "nomsa" && "$template_selection_method" == "standard" ]]; then
            folder_input="json_nomsa_template"
            echo "Using $folder_input"
        else
            folder_input="json_msa_template"
            echo "Using $folder_input"

        fi
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
            $ARRAY_MAP $input_basejson/$folder_input $output_folder $logs_path_datageneration \
            $uniprot_msa $template_selection_method $start)

            echo "Job $jobid finished."
        done

    done
else 
    echo "Skipping AF3 data generation step"
fi

if $RUN_CUSTOM_JSON_GENERATION; then
    echo "Generating custom JSON input files with different MSA and Template settings"
    apptainer exec \
    --bind "$output_datageneration:/mnt/input_data" \
    --bind "$output_customjson:/mnt/output_json" \
    --bind "$src:/mnt/source" \
    "$IMAGE" \
    python /mnt/source/scripts/create_custom_json_reduced.py \
        -i /mnt/input_data \
        -o /mnt/output_json
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
        max_array=1000 
        for start in $(seq $GLOBAL_START $max_array $TOTAL_TASKS); do
            end=$((${start} + ${max_array} - 1))
            [ $end -gt $TOTAL_TASKS ] && end=$TOTAL_TASKS

            array_length=$((${end} - ${start} + 1))
            echo "Submitting array: $start-$end%$CONCURRENT_INFERENCE (length $array_length) for folder: $folder"

            jobid=$(sbatch --parsable --wait \
            --array=1-${array_length}%${CONCURRENT_INFERENCE} \
            scripts/runAF3inference.sh \
            "$folder_path" "$output_inference" "$logs_path_inference" \
            "$ARRAY_MAP" "$NUM_SEEDS" "$NUM_DIFFUSION" "$start")

            echo "Job $jobid for $folder finished."
        done
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
            apptainer exec \
            --bind "$folder_path:/mnt/input_folder" \
            --bind "$TEMPLATE_PATH:/mnt/template_path" \
            --bind "$src:/mnt/source" \
            "$IMAGE" \
            python /mnt/source/scripts/computeDockq.py \
                -i /mnt/input_folder \
                -t /mnt/template_path \
                -n1 D E \
                -m1 D E \
                -n2 C \
                -m2 C

            #python scripts/computeDockq.py -i "$folder_path" -t "$TEMPLATE_PATH" -n1 D E -m1 D E -n2 C -m2 C
            echo "Computation for $folder finished."
        done
        echo "DockQ collected for all folders"
    fi

    echo "Collecting all metrics for AF3 generated structures"
    for folder in $folders_metric_collection; do
        folder_path="$output_inference/$folder"
        echo "Collecting metrics for files in $folder"
        #python scripts/collect_af3metrics_extended.py -i "$folder_path"
        apptainer exec \
        --bind "$folder_path:/mnt/input_folder" \
        --bind "$src:/mnt/source" \
        "$IMAGE" \
        python /mnt/source/scripts/collect_af3metrics_extended.py \
            -i /mnt/input_folder
        echo "Computation for $folder finished."

    done

    echo "Combining metrics in one file"
    apptainer exec \
    --bind "$output_inference:/mnt/output_inference" \
    --bind "$src:/mnt/source" \
    "$IMAGE" \
    python /mnt/source/scripts/combine_metrics_onefile.py \
        -i /mnt/output_inference
    #python combine_metrics_onefile.py -i $output_inference

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
        #folders_nettcrstruc_str="${folders_nettcrstruc[*]}"
        folders_nettcrstruc_str=$(printf '%s\n' "${folders_nettcrstruc[@]}")

        echo "Folders to process: $folders_nettcrstruc_str"
    else
        #folders_nettcrstruc=$FOLDERS_NETTCRSTRUC
        read -r -a folders_nettcrstruc <<< "$FOLDERS_NETTCRSTRUC"
        #folders_nettcrstruc_str="${folders_nettcrstruc[*]}"
        folders_nettcrstruc_str=$(printf '%s\n' "${folders_nettcrstruc[@]}")
        echo "Using folders from config. Folders to process: $folders_nettcrstruc_str"
    fi

    if $PREPARE_NETTCRSTRUCT_INPUT; then
        for folder in $folders_nettcrstruc; do
            folder_path="$output_inference/$folder"
            echo "Submitting nettcrstruct prepare input job for $folder"
            #python scripts/prepare_input_nettcrstruc.py -i $folder_path -o $output_nettcrstruct_datareformatting
            apptainer exec \
            --bind "$folder_path:/mnt/input_folder" \
            --bind "$output_nettcrstruct_datareformatting:/mnt/output_data" \
            --bind "$src:/mnt/source" \
            "$IMAGE" \
            python /mnt/source/scripts/prepare_input_nettcrstruc.py \
                -i /mnt/input_folder \
                -o /mnt/output_data
        done
    else
        echo "Skipping input reformatting for nettcrstruct"
    fi

    num_folders=${#folders_nettcrstruc[@]}
    echo "Submitting netttcrstruct for $num_folders folders"
    jobid=$(sbatch --parsable --wait --array=1-${num_folders}%4 scripts/run_nettcrstruc.sh "$output_nettcrstruct_datareformatting" "$logs_path_nettcrstruct" "$ENSEMBLE" "$suffix_output_inference" "$folders_nettcrstruc_str")
    echo "Job $jobid for nettcrstruct finished."
    #python scripts/collect_nettcrstruct_reranking.py -i $output_nettcrstruct_datareformatting -o $output_inference
    apptainer exec \
    --bind "$output_nettcrstruct_datareformatting:/mnt/input_data" \
    --bind "$output_inference:/mnt/output_inference" \
    --bind "$src:/mnt/source" \
    "$IMAGE" \
    python /mnt/source/scripts/collect_nettcrstruct_reranking.py \
        -i /mnt/input_data \
        -o /mnt/output_inference


else
    echo "Skipping nettcrstruc computation"

fi

echo "Pipeline completed"
