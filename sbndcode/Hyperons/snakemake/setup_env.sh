#!/bin/bash

# ==========================================
# LArSoft Environment Setup
# ==========================================

echo "Sourcing SBND base environment..."
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh


echo "Sourcing local LArSoft products..."
source "$SETUP_LOCAL"

echo "Setting up LArSoft release (${SBNDCODE_VERSION} - ${SBNDCODE_QUALS})..."
setup sbndcode -v $SBNDCODE_VERSION -q $SBNDCODE_QUALS

echo "Injecting local work directory into FHiCL path..."
export FHICL_FILE_PATH="${HYPERON_WORK_DIR}:${FHICL_FILE_PATH}"

echo "Environment and filepath setup complete."