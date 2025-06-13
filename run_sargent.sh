#!/bin/bash
eval "$(conda shell.bash hook)"

set -e

R_ENV_PATH=$1

# Run ENACT pipeline
conda activate $R_ENV_PATH
echo "Running sargent with config: $CONFIG_PATH"
Rscript r_scripts/sargent.R