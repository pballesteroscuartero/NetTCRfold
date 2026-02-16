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

#source configs/af3Benchmark_allfolders_2seeds_50diffusion.cfg
#source configs/schumacher.cfg
#source configs/immrep2025.cfg
#source configs/nettcr.cfg
source configs/boltz2_benchmark.cfg


src=/home/projects2/pbacu/projects/structureTCR/structurePipeline
mhcdb=$src/databases/mhc_sequences
dockq_repo=/home/projects2/pbacu/repositories/DockQ
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
    #combinations=(
    #"paired onquery"
    #"unpaired onquery"
    #)  

    combinations=(
    "paired standard"
    "unpaired standard"
    #"unpaired onquery"
    )  
    
    for combo in "${combinations[@]}"; do
        read -r uniprot_msa template_selection_method <<< "$combo"

        echo "Running AF data generation step with MSA (uniprot in ${uniprot_msa}) and Template with selection method ${template_selection_method}:"  
        output_folder="${output_datageneration}/uniprotOn_${uniprot_msa}_template_${template_selection_method}"
        max_array=1000 
  
        for start in $(seq $GLOBAL_START $max_array $TOTAL_TASKS_DATA); do
            end=$((${start} + ${max_array} - 1))
            [ $end -gt $TOTAL_TASKS_DATA ] && end=$TOTAL_TASKS_DATA

            array_length=$((${end} - ${start} + 1))
            echo "Submitting array: $start-$end%$CONCURRENT (length $array_length) for combo: $combo"

            jobid=$(sbatch --parsable --wait \
            --array=1-${array_length}%${CONCURRENT} \
            scripts/runAF3_dataGeneration.sh \
            $ARRAY_MAP_DATA $input_basejson/"json_msa_template" $output_folder $logs_path_datageneration \
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
    --bind "$INPUT_DB:/mnt/input_db" \
    --bind "$src:/mnt/source" \
    "$IMAGE" \
    python /mnt/source/scripts/create_custom_json_reduced.py \
        -i /mnt/input_data \
        -o /mnt/output_json \
        -d /mnt/input_db/"${INPUT_FILE%.csv}_hla_withid.csv"
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
        for start in $(seq $GLOBAL_START $max_array $TOTAL_TASKS_INFERENCE); do
            end=$((${start} + ${max_array} - 1))
            [ $end -gt $TOTAL_TASKS_INFERENCE ] && end=$TOTAL_TASKS_INFERENCE

            array_length=$((${end} - ${start} + 1))
            echo "Submitting array: $start-$end%$CONCURRENT_INFERENCE (length $array_length) for folder: $folder"

            jobid=$(sbatch --parsable --wait \
            --array=1-${array_length}%${CONCURRENT_INFERENCE} \
            scripts/runAF3inference.sh \
            "$folder_path" "$output_inference" "$logs_path_inference" \
            "$ARRAY_MAP_INFERENCE" "$NUM_SEEDS" "$NUM_DIFFUSION" "$start")

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
            if $COMPUTE_ON_PMHC; then
                apptainer exec \
                    --bind "$folder_path:/mnt/input_folder" \
                    --bind "$TEMPLATE_PATH:/mnt/template_path" \
                    --bind "$src:/mnt/source" \
                    --bind "$dockq_repo:/mnt/dockq_repo" \
                    "$IMAGE" \
                    python /mnt/source/scripts/computeDockq.py \
                        -i /mnt/input_folder \
                        -t /mnt/template_path \
                        -n1 D E \
                        -m1 D E \
                        -n2 C A \
                        -m2 C A \
                        -s "_dockQonpMHC"   
            fi
            if $COMPUTE_ON_PEP; then
                apptainer exec \
                    --bind "$folder_path:/mnt/input_folder" \
                    --bind "$TEMPLATE_PATH:/mnt/template_path" \
                    --bind "$src:/mnt/source" \
                    --bind "$dockq_repo:/mnt/dockq_repo" \
                    "$IMAGE" \
                    python /mnt/source/scripts/computeDockq.py \
                        -i /mnt/input_folder \
                        -t /mnt/template_path \
                        -n1 D E \
                        -m1 D E \
                        -n2 C \
                        -m2 C \
                        -s "_dockQonpeptide"
            fi
            echo "Computation for $folder finished."
            echo "Removing extra pdb files generated"
            find "$folder_path" -mindepth 3 -maxdepth 3 -type f -name "*.pdb" -delete
        done
        echo "DockQ collected for all folders"
    fi

    echo "Collecting all metrics for AF3 generated structures and combining them in one file"

    suffixes=()
    if $COMPUTE_ON_PMHC; then
        suffixes+=("_dockQonpMHC") 
    fi
    if $COMPUTE_ON_PEP; then
        suffixes+=("_dockQonpeptide")  
    fi
    if [ ${#suffixes[@]} -eq 0 ]; then
        suffixes=("")
        echo "DockQ will not be collected, only AF3 metrics will be collected and combined"
    fi

    for folder in $folders_metric_collection; do
        folder_path="$output_inference/$folder"
        echo "Collecting metrics for files in $folder"
        for i in "${!suffixes[@]}"; do
            suffix="${suffixes[$i]}"
        
            apptainer exec \
            --bind "$folder_path:/mnt/input_folder" \
            --bind "$src:/mnt/source" \
            "$IMAGE" \
            python /mnt/source/scripts/collect_af3metrics_extended.py \
                -i /mnt/input_folder \
                -s "$suffix"
        done
        echo "Computation for $folder finished."
    done

    for i in "${!suffixes[@]}"; do
        suffix="${suffixes[$i]}"
        apptainer exec \
                --bind "$output_inference:/mnt/output_inference" \
                --bind "$src:/mnt/source" \
                "$IMAGE" \
                python /mnt/source/scripts/combine_metrics_onefile.py \
                    -i /mnt/output_inference \
                    -s "$suffix"
    done

else
    echo "Skipping metrics collection"

    
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
        for folder in $folders_nettcrstruc_str; do
            folder_path="$output_inference/$folder"
            mkdir -p $output_nettcrstruct_datareformatting
            echo "Submitting nettcrstruct prepare input job for $folder"
            #python scripts/prepare_input_nettcrstruc.py -i $folder_path -o $output_nettcrstruct_datareformatting
            apptainer exec \
            --bind "$folder_path:/mnt/input_folder" \
            --bind "$output_nettcrstruct_datareformatting:/mnt/output_data" \
            --bind "$src:/mnt/source" \
            "$IMAGE" \
            python /mnt/source/scripts/prepare_input_nettcrstruc.py \
                -i /mnt/input_folder \
                -o /mnt/output_data \
                -n $folder
        done
    else
        echo "Skipping input reformatting for nettcrstruct"
    fi

    num_folders=${#folders_nettcrstruc[@]}
    echo "Submitting netttcrstruct for $num_folders folders"
    jobid=$(sbatch --parsable --wait --array=1-${num_folders}%4 scripts/run_nettcrstruc.sh "$output_nettcrstruct_datareformatting" "$logs_path_nettcrstruct" "$ENSEMBLE" "$suffix_output_inference" "$folders_nettcrstruc_str")
    echo "sbatch exit code: $?"
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
    echo "Removing extra cif files generated"
    find "$output_nettcrstruct_datareformatting" -mindepth 3 -maxdepth 3 -type f -name "*.cif" -delete


else
    echo "Skipping nettcrstruc computation"

fi

echo "Pipeline completed"
