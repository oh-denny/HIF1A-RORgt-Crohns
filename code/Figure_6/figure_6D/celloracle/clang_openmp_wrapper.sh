#!/usr/bin/env bash
set -euo pipefail

CONDA_PREFIX="${CONDA_PREFIX:-/Users/denny/miniconda3/envs/celloracle-hif1a}"

args=()
for arg in "$@"; do
  if [[ "$arg" == "-fopenmp" ]]; then
    args+=("-Xpreprocessor" "-fopenmp")
  else
    args+=("$arg")
  fi
done

exec /usr/bin/clang \
  -I"${CONDA_PREFIX}/include" \
  -L"${CONDA_PREFIX}/lib" \
  "${args[@]}" \
  -lomp
