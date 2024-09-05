import numpy as np
import scipy.stats
from .lsalib_stats import tied_rank
def percentileNormalize(tseries):
    ranks = tied_rank(tseries)
    nt = np.ma.masked_invalid(scipy.stats.norm.ppf(ranks / (len(ranks) - np.sum(ranks.mask) + 1)), copy=True)
    nt = nt.filled(fill_value=0)
    return nt

def percentileZNormalize(tseries):
    ranks = tied_rank(tseries)
    nt = np.ma.masked_invalid(scipy.stats.norm.ppf(ranks / (len(ranks) - np.sum(ranks.mask) + 1)), copy=True)
    try:
        zt = (nt - np.ma.mean(nt, axis=0)) / np.ma.std(nt)
    except FloatingPointError:
        zt = nt - np.ma.mean(nt, axis=0)
    zt = zt.filled(fill_value=0)
    return zt

def robustZNormalize(tseries):
    ranks = tied_rank(tseries)
    nt = np.ma.masked_invalid(scipy.stats.norm.ppf(ranks / (len(ranks) - np.sum(ranks.mask) + 1)), copy=True)
    mad_sd = 1.4826 * np.ma.median(np.ma.abs(nt - np.ma.median(nt)))
    range_sd = (np.ma.max(nt) - np.ma.min(nt)) / 4
    sd_est = range_sd if mad_sd == 0 else mad_sd
    try:
        zt = (nt - np.ma.median(nt)) / sd_est
    except FloatingPointError:
        zt = nt - np.ma.median(nt)
    zt = zt.filled(fill_value=0)
    return zt

def noZeroNormalize(tseries):
    nt = np.ma.masked_equal(tseries, 0)
    if type(nt.mask) == np.bool_:
        nt.mask = [nt.mask] * nt.shape[0]
    ranks = tied_rank(nt)
    nt = np.ma.masked_invalid(scipy.stats.norm.ppf(ranks / (len(ranks) - np.sum(ranks.mask) + 1)), copy=True)
    try:
        zt = (nt - np.ma.mean(nt, axis=0)) * (1 / np.ma.std(nt, axis=0))
    except FloatingPointError:
        zt = nt - np.ma.mean(nt, axis=0)
    zt = np.ma.masked_invalid(zt)
    zt = zt.filled(fill_value=0)
    return zt

def noneNormalize(tseries):
    nt = tseries.filled(fill_value=0)
    return nt

# Add other normalization functions here