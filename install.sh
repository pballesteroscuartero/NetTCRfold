#!/bin/bash
# Condensed installation steps from README.md's "Installation" section.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "== Step 2: AF3 weights =="
echo "Manual step — this cannot be scripted. Request access to the AF3 weights"
echo "and place them inside $SCRIPT_DIR/alphafold3_tcrpmhc (see README.md)."
echo

echo "== Step 3: AF3 pipeline image =="
wget -P apptainer \
    https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/image/

echo "== Step 4: Chemical components =="
wget -P alphafold3_tcrpmhc/src/alphafold3/constants/converters \
    https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/chemicalComponents/ccd.pickle \
    https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/chemicalComponents/chemical_component_sets.pickle

echo "== Step 5: Curated databases =="
wget -r -np -nH --cut-dirs=3 -R "index.html*" -P alphafold3_tcrpmhc/ \
    https://services.healthtech.dtu.dk/suppl/immunology/NetTCRaFold-1.0/tcrpmhc_databases/

echo "== Step 6: DockQ =="

git clone git@github.com:wallnerlab/DockQ.git
cd DockQ
git checkout 3735c16

echo "== Step 7: Conda environment =="
conda env create -f pipeline/NetTCRfold.yml
conda activate NetTCRfold
conda run -n NetTCRfold pip install -e "$SCRIPT_DIR/pipeline"

echo
echo "Installation finished"
