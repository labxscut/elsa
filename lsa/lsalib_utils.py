import numpy as np
import scipy.interpolate
import subprocess
import os

def ma_average(ts, axis=0):
    ns = np.ma.mean(ts, axis=0)
    if type(ns.mask) == np.bool_:
        ns.mask = [ns.mask] * ns.shape[axis]
    return ns

def ma_median(ts, axis=0):
    ns = np.ma.median(ts, axis=axis)
    if type(ns.mask) == np.bool_:
        ns.mask = [ns.mask] * ns.shape[axis]
    return ns

def sample_wr(population, k):
    n = len(population)
    _random, _int = np.random.random, int
    result = np.array([np.nan] * k)
    for i in range(k):
        j = _int(_random() * n)
        if type(population) == np.ma.MaskedArray:
            if population.mask[j]:
                result[i] = np.nan
            else:
                result[i] = population[j]
        else:
            result[i] = population[j]
    if type(population) == np.ma.MaskedArray:
        result = np.ma.masked_invalid(result)
    return result

def fillMissing(tseries, method):
    if method == 'none':
        return tseries
    
    valid = ~np.isnan(tseries)
    indices = np.arange(len(tseries))
    
    if method == 'zero':
        tseries[~valid] = 0
    elif method in ['linear', 'quadratic', 'cubic']:
        f = scipy.interpolate.interp1d(indices[valid], tseries[valid], kind=method, bounds_error=False, fill_value='extrapolate')
        tseries = f(indices)
    elif method == 'nearest':
        f = scipy.interpolate.interp1d(indices[valid], tseries[valid], kind='nearest', bounds_error=False, fill_value='extrapolate')
        tseries = f(indices)
    
    return tseries

def safeCmd(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()
    except subprocess.CalledProcessError as e:
        return f"Error: {e}"

def float_equal(a, b, tol=1e-6):
    return abs(a - b) <= tol

# Add any other utility functions here that were in the original lsalib.py
# but not included in other new modules