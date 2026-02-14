#!/bin/bash
#SBATCH --job-name=fragfold_watcher
#SBATCH -n 4
#SBATCH --output=logs/fragfold_%A_%a.out
#SBATCH --error=logs/fragfold_%A_%a.err
#SBATCH --partition=pi_keating

MAX_JOBS=2
PARAM_FILE="./params.yaml"

source activate fragfold3
fragfold3 --input_params "$PARAM_FILE" --colabfold_scheduler slurm --colabfold_max_jobs_allowed "$MAX_JOBS"





