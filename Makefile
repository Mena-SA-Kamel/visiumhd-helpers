ENV_DIR := /home/oneai/envs/

R_ENV_NAME := enact_r_env

CONFIG_PATH ?= /home/oneai/enact-pipeline/config/configs.yaml

R_ENV_PATH := $(ENV_DIR)$(R_ENV_NAME)

run_sargent:
	bash run_sargent.sh $(R_ENV_PATH)

setup_r_env:
	bash setup_r_env.sh $(R_ENV_PATH)
