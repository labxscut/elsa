"""
Compatibility layer for compcore module.
This module provides a bridge between the old SWIG interface and new PyBind11 bindings.
It re-exports all the necessary symbols from _compcore while maintaining backward compatibility.

Classes:
    LSA_Data: Container for Local Similarity Analysis input data
    LSA_Result: Container for Local Similarity Analysis results
    LLA_Data: Container for Local Liquid Association input data
    LLA_Result: Container for Local Liquid Association results

Functions:
    DP_lsa: Compute Local Similarity Analysis using dynamic programming
    DP_lla: Compute Local Liquid Association using dynamic programming
    calc_LA: Calculate static Liquid Association
    test: Simple test function
"""

# Import all symbols from the compiled _compcore module
from ._compcore import (
    LSA_Data,      # Container for LSA input data
    LSA_Result,    # Container for LSA results
    DP_lsa,        # Main LSA computation function
    LLA_Data,      # Container for LLA input data
    LLA_Result,    # Container for LLA results
    DP_lla,        # Main LLA computation function
    calc_LA,       # Static LA calculation
    test           # Test function
)

# Make these symbols available at module level for backward compatibility
__all__ = [
    'LSA_Data',
    'LSA_Result',
    'DP_lsa',
    'LLA_Data',
    'LLA_Result',
    'DP_lla',
    'calc_LA',
    'test'
] 