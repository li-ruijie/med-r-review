#!/usr/bin/env bash
# Launches 07_simulation.R detached from the terminal (nohup).
#
# Usage (from the directory that holds code/): bash code/run_simulation.sh pilot|check|full [method ...]
#
# Naming methods re-runs only those methods in the finished scenarios (see
# 07_simulation.R).
#
# BLAS is pinned to one thread per process before R starts, since the
# simulation already runs one worker per hardware thread.

set -euo pipefail

mode="${1:-pilot}"
shift || true
cd "$(dirname "$0")"
out_dir="../results/simulation/${mode}"
mkdir -p "${out_dir}"

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

nohup Rscript 07_simulation.R "${mode}" "$@" >> "${out_dir}/console.log" 2>&1 < /dev/null &
echo "started ${mode} (pid $!), log in ${out_dir}"
