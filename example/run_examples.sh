#!/usr/bin/env bash
# Examples for Dataset_collection.py and dataset_script.py
#
# Prerequisites:
#   pip install rdkit pandas requests
#   Download ChEMBL SQLite DB (optional, for local use):
#     https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/
#
# All examples below assume you run them from the repo root directory.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
SMI_STD="$SCRIPT_DIR/chembl_33_step4.smi.tar.gz"

# ---------------------------------------------------------------------------
# Dataset_collection.py  — single target, online ChEMBL API
# ---------------------------------------------------------------------------

# Main example: CHEMBL301 (Carbonic anhydrase II), regression, with standardized SMILES
python Dataset_collection.py \
    -i CHEMBL301 \
    --t reg \
    -v plog_value \
    -s "$SMI_STD" \
    --activity Ki IC50 EC50 Kd \
    --active_value 6 --active_op ">=" \
    --inactive_value 5 --inactive_op "<="

# Classification dataset: CHEMBL6030 (enzyme), Ki/IC50, plog_value threshold
python "$REPO_DIR/Dataset_collection.py" \
    -i CHEMBL6030 \
    --t class \
    --activity Ki IC50 \
    -s "$SMI_STD" \
    -a allosteric \
    --active_value 7 --active_op ">=" \
    --inactive_value 5 --inactive_op "<" \
    -o output/CHEMBL6030/CHEMBL6030

# Using a local ChEMBL SQLite database instead of the API
# python "$REPO_DIR/Dataset_collection.py" \
#     -i CHEMBL6030 \
#     --t class \
#     --activity Ki IC50 \
#     -d /path/to/chembl_33.db \
#     -s "$SMI_STD" \
#     -o output/CHEMBL6030/CHEMBL6030

# Show the built-in activity-type classification tree and exit
python "$REPO_DIR/Dataset_collection.py" \
    -i CHEMBL6030 --t class --activity Ki \
    --show_tree_act

# ---------------------------------------------------------------------------
# dataset_script.py  — batch processing from a CSV file
# ---------------------------------------------------------------------------

# Using target_class column (applies class-specific thresholds from Bosc et al.)
python "$REPO_DIR/dataset_script.py" \
    -i "$SCRIPT_DIR/chembl29_targets_test.csv" \
    -t class \
    -s "$SMI_STD"

# Using custom per-target cutoffs (cutoff_act / cutoff_inact columns)
python "$REPO_DIR/dataset_script.py" \
    -i "$SCRIPT_DIR/targets_custom_cutoffs.csv" \
    -t class \
    -s "$SMI_STD"

# Restricting activity types via the category column
python "$REPO_DIR/dataset_script.py" \
    -i "$SCRIPT_DIR/targets_category.csv" \
    -t class \
    -s "$SMI_STD"

# Regression mode with 4 parallel workers and standardized SMILES
python "$REPO_DIR/dataset_script.py" \
    -i "$SCRIPT_DIR/chembl29_targets_test.csv" \
    -t reg \
    -n 4 \
    -s "$SMI_STD"

# With a local ChEMBL SQLite database
# python "$REPO_DIR/dataset_script.py" \
#     -i "$SCRIPT_DIR/chembl29_targets_test.csv" \
#     -t class \
#     -d /path/to/chembl_33.db \
#     -s "$SMI_STD" \
#     -n 4
