#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -J sampling
#SBATCH --array=0-51
# One (variant, scenario) per array task: 4 variants x 13 scenarios.
# Usage: sbatch paper_figures/sample_array.sh [optFrac]   (default 0.99999, the optimal face)
# Needs RAVEN_ROOT, GUROBI_HOME, POLYTOPESAMPLER_ROOT, MATLAB_PREFDIR and YEAST_GEM_ROOT in the environment.
set -euo pipefail
VARIANTS=(v910_new v910_old v902_new v902_old)
SCENARIOS=(sjoberg bakker_wt bakker_adh3 bakker_adh3adh4 boles_wt boles_mae1 flux_Celton2012 flux_Jouhten2008 flux_Nissen01 flux_Nissen02 flux_Nissen03 flux_Nissen04 flux_Wasylenko2014)
V=${VARIANTS[$((SLURM_ARRAY_TASK_ID / 13))]}
SC=${SCENARIOS[$((SLURM_ARRAY_TASK_ID % 13))]}
OPT=${1:-0.99999}
cd "$(dirname "$0")/../r1_major5_degeneracy"
matlab -batch "runSampling('$V', '$SC', 1000, 200, 1, $OPT)"
