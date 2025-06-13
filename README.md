# visiumhd-helpers
R helper functions for Visium HD analysis. This can be helpful to run R-based packages that are used by packages such as ENACT (https://github.com/Sanofi-Public/enact-pipeline.git).

Running Sargent on ENACT bin-to-cell assignment results and get compatible cell annotation results.

Install
1. Update Makefile with:
ENV_DIR -> direcotry to save the R conda environment
R_ENV_NAME -> R environment name
CONFIG_PATH => path to enact config file

2. Run make setup_r_env

Run Sargent on data.
Run make run_sargent
