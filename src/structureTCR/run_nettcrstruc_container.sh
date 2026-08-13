#! /bin/bash
# static variables (expect alphafold3 resources to be in directory as this file)
export NTS_RESOURCES_DIR="/home/projects2/pbacu/repositories/NetTCR-struc/nettcrstruc"
export NTS_SRC=${NTS_RESOURCES_DIR}
#export NTS_IMAGE=${NTS_RESOURCES_DIR}/image/nettcrstruct_python310_cuda126_forGefion_withSymlink.sif
export NTS_IMAGE=/home/projects2/pbacu/repositories/containers/NetTCR-struc/image/nettcrstruct_python310_cuda117_withSymlink_TCRpepRanking.sif

#Input variables
export NTS_INPUTFOLDER="${1%/}"
export NTS_FEATURESSAVE=$2
export NTS_ENSEMBLE_GVP_IF1=$3
export NTS_ENSEMBLE_GVP=$4

NTS_FEATURES_LOAD_GVP_IF1=$NTS_FEATURESSAVE/gvp_if1_embeddings
NTS_FEATURES_LOAD_GVP=$NTS_FEATURESSAVE/gvp

mkdir -p "$NTS_FEATURESSAVE"
mkdir -p "$NTS_FEATURES_LOAD_GVP_IF1"
mkdir -p "$NTS_FEATURES_LOAD_GVP"
# print paths 
echo NetTCRStruct resource: $NTS_RESOURCES_DIR
echo NetTCRStruct source directory $NTS_SRC
echo NetTCRStruct image: $NTS_IMAGE
echo NetTCRStruct input: $NTS_INPUTFOLDER
echo NetTCRStruct features save dir: $NTS_FEATURESSAVE
echo NetTCRStruct ensemble gvp: $NTS_ENSEMBLE_GVP
echo NetTCRStruct ensemble gvp if1: $NTS_ENSEMBLE_GVP_IF1

#REAL_PARENT=$(realpath "${NTS_INPUTFOLDER}/..")

echo "Creating features"
apptainer exec \
     --nv \
     --bind $NTS_INPUTFOLDER:/mnt/nts_input \
     --bind $NTS_FEATURESSAVE:/mnt/nts_features \
     --bind $NTS_SRC:/mnt/nts_source/ \
     $NTS_IMAGE \
     python /mnt/nts_source/scripts/create_geometric_features.py \
     -i /mnt/nts_input \
     -o /mnt/nts_features \
     -n 2 \
     -d cuda \
     --chain_names D E C A
     
echo "Reranking with gnn gvp if1"
apptainer exec \
     --nv \
     --bind $NTS_INPUTFOLDER:/mnt/nts_input \
     --bind $NTS_FEATURES_LOAD_GVP_IF1:/mnt/nts_features_gvp_if1 \
     --bind $NTS_SRC:/mnt/nts_source \
     $NTS_IMAGE \
     python /mnt/nts_source/scripts/rerank_docking_poses_withTCRpep.py \
     input_dir=/mnt/nts_input \
     processed_dir=/mnt/nts_features_gvp_if1 \
     name=ensemble_gvp_if1 \
     ensemble=$NTS_ENSEMBLE_GVP_IF1 \
     chain_names=[D,E,C,A]

echo "Reranking with gnn gvp"
apptainer exec \
     --nv \
     --env PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:128 \
     --bind $NTS_INPUTFOLDER:/mnt/nts_input \
     --bind $NTS_FEATURES_LOAD_GVP:/mnt/nts_features_gvp \
     --bind $NTS_SRC:/mnt/nts_source \
     $NTS_IMAGE \
     python /mnt/nts_source/scripts/rerank_docking_poses_withTCRpep.py \
     input_dir=/mnt/nts_input \
     processed_dir=/mnt/nts_features_gvp \
     name=ensemble_gvp \
     ensemble=$NTS_ENSEMBLE_GVP \
     chain_names=[D,E,C,A] \
    +ensemble.model_0.model.node_in_dim=[30,3] \
    +ensemble.model_1.model.node_in_dim=[30,3] \
    +ensemble.model_2.model.node_in_dim=[30,3] \
    +ensemble.model_3.model.node_in_dim=[30,3] \
    +ensemble.model_4.model.node_in_dim=[30,3] \
    +ensemble.model_5.model.node_in_dim=[30,3] \
    +ensemble.model_6.model.node_in_dim=[30,3] \
    +ensemble.model_7.model.node_in_dim=[30,3] \
    +ensemble.model_8.model.node_in_dim=[30,3] \
    +ensemble.model_9.model.node_in_dim=[30,3]

#     --bind ${REAL_PARENT}:${REAL_PARENT} \
