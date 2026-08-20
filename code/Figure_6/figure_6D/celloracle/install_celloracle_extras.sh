#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-celloracle-hif1a}"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WRAPPER="${PROJECT_DIR}/code/celloracle/clang_openmp_wrapper.sh"

chmod +x "${WRAPPER}"

CC="${WRAPPER}" CXX="${WRAPPER}" \
  conda run -n "${ENV_NAME}" \
  python -m pip install --no-build-isolation "velocyto>=0.17"

conda run -n "${ENV_NAME}" \
  python -m pip install --no-deps "celloracle==0.20.0"
