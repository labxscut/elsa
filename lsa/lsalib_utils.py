import numpy as np
import scipy.interpolate

def ma_average(ts, axis=0):
    ns = np.ma.mean(ts, axis=0)
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

# Add other utility functions here