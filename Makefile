ENV_DIR := <conda_env_directory>

R_ENV_NAME := r_env

CONFIG_PATH ?= <path_to_config_file>

R_ENV_PATH := $(ENV_DIR)$(R_ENV_NAME)


setup_r_env:
	bash setup_r_env.sh $(R_ENV_PATH)

run_sargent:
	CONFIG_PATH=$(CONFIG_PATH) bash run_sargent.sh $(R_ENV_PATH)
