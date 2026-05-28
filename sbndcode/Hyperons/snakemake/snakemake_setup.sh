#!/bin/bash

# ==========================================
# Background Paths
# ==========================================
export HYPERON_BG_INPUT="/data/tjones9/SBND/2026_production" 
export HYPERON_BG_OUTPUT="/data/tjones9/sbnd/hyperons/hyperon_analyser_output/2026_production"

# ==========================================
# Signal Paths
# ==========================================
export HYPERON_SIG_INPUT="/data/sbnd/hyperons_new/*batch*"
# Mirrors the batch structure into your main output directory
export HYPERON_SIG_OUTPUT="/data/tjones9/sbnd/hyperons/hyperon_analyser_output/filtered_hyperons"

# ==========================================
# Global Project Config
# ==========================================
export HYPERON_WORK_DIR="$HOME/SBND/larsoft_v10_14_02_03/srcs/sbndcode/sbndcode/Hyperons"
export SETUP_LOCAL=${1:-"$HOME/SBND/larsoft_v10_14_02_03/localProducts_larsoft_v10_14_02_03_prof_e26/setup"}
export LOCAL_PRODS="$HOME/SBND/larsoft_v10_14_02_03/localProducts_larsoft_v10_14_02_03_prof_e26"

export CONTAINER_IMAGE="/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest"
export APPTAINER_BIN="/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer"

export SBNDCODE_VERSION="v10_14_02_03"
export SBNDCODE_QUALS="e26:prof"