#!/usr/bin/env bash
set -euo pipefail

workspace="/workspace/elsa"
cd "${workspace}"

gpu_so_count="$(find lsa -maxdepth 1 -name '_compcore_gpu*.so' | wc -l | tr -d ' ')"
if [[ "${gpu_so_count}" == "0" ]]; then
  echo "[elsa-entrypoint] GPU extension not found in mounted workspace; building it now." >&2
  bash ./build_gpu.sh
fi

exec "$@"
