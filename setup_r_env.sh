#!/bin/bash
eval "$(conda shell.bash hook)"

set -e

R_ENV_PATH=$1

# Create R environment
if ! conda info --envs | grep -q "$R_ENV_PATH"; then
    echo "Environment $R_ENV_PATH does not exist. Creating..."
    conda config --set ssl_verify False
    conda create --prefix $R_ENV_PATH r-essentials r-irkernel r-rcpparmadillo r-devtools
    conda activate $R_ENV_PATH
    sudo apt update
    sudo apt install build-essential
    conda install conda-forge::r-ggdist
    conda install -c conda-forge r-seurat -y
    echo "Installing sargent"
    Rscript -e 'remotes::install_git("https://github.com/nourin-nn/sargent.git", build_vignettes = FALSE)'
fi