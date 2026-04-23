"""Backend selector for the ``lsa.compcore`` compatibility layer.

The repository currently ships with the CPU extension built as ``lsa._compcore``.
This selector keeps that path as the default and will optionally prefer a future
GPU extension named ``lsa._compcore_gpu`` when requested via the environment.

Environment variables:
    ELSA_COMPCORE_BACKEND:
        - ``auto`` (default): try GPU first, then fall back to CPU
        - ``gpu``: require the GPU backend if present, otherwise fall back to CPU
        - ``cpu``: force the CPU backend
"""

from __future__ import annotations

import os
from importlib import import_module
from types import ModuleType


def _resolve_backend_name() -> str:
    backend = os.getenv("ELSA_COMPCORE_BACKEND", "auto").strip().lower()
    if backend in {"cpu", "gpu", "auto"}:
        return backend
    return "auto"


def _import_candidate(module_name: str) -> ModuleType | None:
    try:
        return import_module(module_name, package=__package__)
    except ImportError:
        return None


def _load_backend() -> tuple[str, ModuleType]:
    requested_backend = _resolve_backend_name()
    candidates: list[tuple[str, str]]

    if requested_backend == "cpu":
        candidates = [("cpu", "._compcore")]
    elif requested_backend == "gpu":
        candidates = [("gpu", "._compcore_gpu"), ("cpu", "._compcore")]
    else:
        candidates = [("gpu", "._compcore_gpu"), ("cpu", "._compcore")]

    for backend_name, module_name in candidates:
        module = _import_candidate(module_name)
        if module is not None:
            return backend_name, module

    raise ImportError("Unable to import any compcore backend")


BACKEND_NAME, _backend = _load_backend()
USING_GPU = BACKEND_NAME == "gpu"

LSA_Data = _backend.LSA_Data
LSA_Result = _backend.LSA_Result
DP_lsa = _backend.DP_lsa
LLA_Data = _backend.LLA_Data
LLA_Result = _backend.LLA_Result
DP_lla = _backend.DP_lla
calc_LA = _backend.calc_LA
test = _backend.test
