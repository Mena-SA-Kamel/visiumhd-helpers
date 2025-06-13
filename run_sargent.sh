#!/bin/bash
eval "$(conda shell.bash hook)"

set -e

R_ENV_PATH=$1

# Run ENACT pipeline
conda activate $R_ENV_PATH
echo "run sargent"
Rscript src/enact/Sargent.R