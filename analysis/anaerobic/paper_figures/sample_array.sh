#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -J sampling
#SBATCH --array=0-147
# One RHMC chain per array task. Tasks 0-51: seed 1 of 4 models x 13 scenarios.
# Tasks 52-147: seeds 2-4 of 4 models x 8 flux scenarios, for the convergence
# diagnostics (chain_diagnostics.py).
# Usage, from analysis/anaerobic: sbatch [-A account -t 08:00:00] paper_figures/sample_array.sh [optFrac]
#   optFrac 0.99999 (default) samples the optimal face, 1.01 a 1% band.
# The environment (MATLAB and Gurobi modules, RAVEN_ROOT, GUROBI_HOME,
# POLYTOPESAMPLER_ROOT, MATLAB_PREFDIR, YEAST_GEM_ROOT) is that of the caller.
set -euo pipefail
VARIANTS=(vdev_new vdev_old v902_new v902_old)
SCENARIOS=(sjoberg bakker_wt bakker_adh3 bakker_adh3adh4 boles_wt boles_mae1 flux_Celton2012 flux_Jouhten2008 flux_Nissen01 flux_Nissen02 flux_Nissen03 flux_Nissen04 flux_Wasylenko2014)
FLUX=(sjoberg flux_Celton2012 flux_Jouhten2008 flux_Nissen01 flux_Nissen02 flux_Nissen03 flux_Nissen04 flux_Wasylenko2014)
T=$SLURM_ARRAY_TASK_ID
if [ "$T" -lt 52 ]; then
    V=${VARIANTS[$((T / 13))]}; SC=${SCENARIOS[$((T % 13))]}; SEED=1
else
    U=$((T - 52))
    V=${VARIANTS[$((U / 24))]}; SC=${FLUX[$(((U % 24) / 3))]}; SEED=$((2 + U % 3))
fi
OPT=${1:-0.99999}
cd "$SLURM_SUBMIT_DIR/r1_major5_degeneracy"
echo "Threads 1" > gurobi.env
matlab -batch "addpath('../common'); runSampling('$V', '$SC', 1000, 200, $SEED, $OPT)"
