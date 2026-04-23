#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python_include="$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["include"])')"
pybind_include="$(python3 -c 'import pybind11; print(pybind11.get_include())')"
extension_suffix="$(python3-config --extension-suffix)"
cuda_home="${CUDA_HOME:-/usr/local/cuda}"

output_dir="${root_dir}/lsa"
gpu_object="${output_dir}/compcore_gpu.o"
gpu_module="${output_dir}/_compcore_gpu${extension_suffix}"

rm -f "${gpu_object}" "${gpu_module}"

nvcc -std=c++14 -Xcompiler -fPIC \
  -I"${root_dir}/lsa" \
  -I"${python_include}" \
  -I"${pybind_include}" \
  -c "${output_dir}/compcore_gpu.cu" \
  -o "${gpu_object}"

g++ -std=c++14 -shared -fPIC \
  -I"${root_dir}/lsa" \
  -I"${python_include}" \
  -I"${pybind_include}" \
  "${output_dir}/compcore_gpu_bindings.cpp" \
  "${gpu_object}" \
  -L"${cuda_home}/lib64" \
  -lcudart \
  $(python3-config --ldflags) \
  -o "${gpu_module}"

rm -f "${gpu_object}"
