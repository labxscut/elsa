"""Compatibility layer for compcore module.

This module keeps the historical import path stable while delegating symbol
resolution to a backend selector. The default backend is the existing CPU
extension. If a GPU extension is added later, it can be selected through the
``ELSA_COMPCORE_BACKEND`` environment variable without changing callers.
"""

from .compcore_backend import (  # noqa: F401
    BACKEND_NAME,
    USING_GPU,
    LSA_Data,
    LSA_Result,
    DP_lsa,
    LLA_Data,
    LLA_Result,
    DP_lla,
    calc_LA,
    test,
)

__all__ = [
    'BACKEND_NAME',
    'USING_GPU',
    'LSA_Data',
    'LSA_Result',
    'DP_lsa',
    'LLA_Data',
    'LLA_Result',
    'DP_lla',
    'calc_LA',
    'test',
]