import numpy as np
import scipy.stats
from scipy import interpolate
import sys
from .lsalib_core import Q_lam_max, Q_lam_step, my_decimal, Rmax_min, Rmax_max, kcut_min, pipi, pipi_inv

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

def float_equal(a, b, tol=1e-6):
    return abs(a - b) <= tol

def tied_rank(values):
    assert type(values) == np.ma.MaskedArray
    V = np.ma.asarray(values)
    nans = (np.nonzero(V.mask)[0]).tolist()
    V = V[~V.mask]
    
    v_num = {}
    v_cum = {}
    for v in V:
        v_num[v] = v_num.get(v, 0) + 1
    suvs = list(v_num.keys())
    suvs.sort()
    c = 0
    for v in suvs:
        c += v_num[v]
        v_cum[v] = c
    
    sV = np.array([(2*v_cum[V[i]]-v_num[V[i]]+1)/2 for i in range(len(V))], dtype='float')
    
    for idx in nans:
        sV = np.insert(sV, idx, np.nan)
    sV = np.ma.masked_invalid(sV, copy=True)
    
    return sV

def simpleMedian(tseries):
    """ simple median

    Args:
      tseries(2d np.mp.array):  one time series with replicates, each row is a replicate

    Returns:
      1d np.mp.array: one row with replicates summarized by median

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
    """
    
    Xf = ma_median(tseries, axis=0)
    return Xf

def madMedian(tseries):
    """	MAD weighted averaging 

    Args:
      tseries(2d np.ma.array):  one time series with replicates, each row is a replicate

    Returns:
      1d np.ma.array: one row with replicates summarized by MAD weighted median

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
    """
    Xf = tseries
    mad = ma_median(np.ma.abs(Xf - ma_median(Xf, axis=0)), axis=0)
    if np.any(mad.mask) or (np.ma.sum(mad==0))>0:
        return simpleMedian(tseries)                  #mad = 0, fall back to simpleMedian
    Xf = ma_median(Xf, axis=0)*(1/mad)*(1/np.ma.sum(1/mad))*(1/mad)                   #mad-weighted sample
    return Xf

def simpleAverage(tseries):
    """ simple averaging 

    Args:
      tseries(np.ma.array):  one 2d time series (masked array) with replicates, each row is a replicate

    Returns:
      (1d np.ma.array) one row with replicates averaged

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
    """
    Xf = ma_average(tseries, axis=0)
    return Xf

def sdAverage(tseries):
    """	SD weighted averaging 

    Args:
      tseries(np.ma.array):  one 2d time series (masked array) with replicates, each row is a replicate

    Returns:
      (1d np.ma.array): one row with replicates SD weighted averaged

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
    """
    try:
        sd = np.ma.std(tseries, axis=0, ddof=1)
    except FloatingPointError:
        return simpleAverage(tseries)                       #sd = 0, fall back to simpleAverage
    if np.any(sd.mask) or (np.ma.sum(sd==0))>0:
        return simpleAverage(tseries)                       #sd = 0, fall back to simpleAverage
    Xf = ma_average(tseries, axis=0)*(1/sd)*(1/np.ma.sum(1/sd))*(1/sd)   #sd-weighted sample
    return Xf

def storeyQvalue(pvalues, lam=np.arange(0, Q_lam_max, Q_lam_step), method='smoother', robust=False, smooth_df=3):
    """ do Q-value calculation

    Args:
      pvalues(np.array):  a set of p-values
      lam(np.array):  tentative lambda data
      method(str):  calculating method, currently only support 'smoother'
      robust(bool): use robust static or not, default not
      smooth_df(int): order of spline function

    Returns:
      qvalues(np.array): a set of qvalues
    """
    
    try:
        mpvalues = np.ma.masked_invalid(pvalues, copy=True)
        rpvalues = mpvalues[~mpvalues.mask]
        p_num = len(pvalues)
        rp_num = len(rpvalues)

        if rp_num <= 1:
            return np.array([np.nan] * p_num, dtype='float')

        rp_max = np.max(rpvalues)
        rp_lam = lam[lam < rp_max]

        if len(rp_lam) <= 1:
            return np.array([np.nan if np.isnan(pvalues[i]) else 0 for i in range(p_num)], dtype='float')

        pi_set = np.zeros(len(rp_lam), dtype='float')
        for i in range(len(rp_lam)):
            pi_set[i] = np.mean(rpvalues >= rp_lam[i]) / (1 - rp_lam[i])

        if method == 'smoother':
            spline_fit = interpolate.interp1d(rp_lam, pi_set, kind=smooth_df)
            pi_0 = spline_fit(np.max(rp_lam))
            pi_0 = np.max([np.min([np.min(pi_0), 1]), 0])
            if pi_0 == 0:
                method = 'bootstrap'

        if method == 'bootstrap':
            pi_min = np.min(pi_set)
            mse = np.zeros((100, len(rp_lam)), dtype='float')
            pi_set_boot = np.zeros((100, len(rp_lam)), dtype='float')
            for j in range(100):
                p_boot = np.random.choice(rpvalues, rp_num, replace=True)
                for i in range(len(rp_lam)):
                    pi_set_boot[j][i] = np.mean(p_boot >= rp_lam[i]) / (1 - rp_lam[i])
                mse[j] = (pi_set_boot[j] - pi_min)**2
            min_mse_j = np.argmin(mse)
            pi_0 = np.min(pi_set_boot[min_mse_j])
            pi_0 = np.max([np.min([np.min(pi_0), 1]), 0])
            if pi_0 == 0:
                pi_0 = Q_lam_step

        rp_argsort = np.argsort(rpvalues)
        rp_ranks = tied_rank(rpvalues)
        if robust:
            rqvalues = pi_0 * rp_num * rpvalues * (1 / (rp_ranks * (1 - np.power((1 - rpvalues), rp_num))))
        else:
            rqvalues = pi_0 * rp_num * rpvalues * (1 / rp_ranks)
        
        rqvalues[rp_argsort[rp_num-1]] = np.min([rqvalues[rp_argsort[rp_num-1]], 1])
        for i in reversed(range(rp_num-1)):
            rqvalues[rp_argsort[i]] = np.min([rqvalues[rp_argsort[i]], rqvalues[rp_argsort[i+1]], 1])

        qvalues = np.array([np.nan] * p_num)
        j = 0
        for i in range(p_num):
            if not mpvalues.mask[i]:
                qvalues[i] = rqvalues[j]
                j += 1

    except:
        print("caution: q-value estimation error", file=sys.stderr)
        print("from scipy: unusable pvalues -> ", rpvalues, file=sys.stderr)
        qvalues = np.array([np.nan] * p_num, dtype='float')

    return qvalues

def scipy_spearmanr(Xz, Yz):
    try:
        return scipy.stats.spearmanr(Xz, Yz)
    except:
        return (np.nan, np.nan)

def calc_spearmanr(Xz, Yz):
    mask = np.logical_or(Xz.mask, Yz.mask)
    Xz.mask = mask
    Yz.mask = mask
    (SCC, P_SCC) = scipy_spearmanr(Xz.compressed(), Yz.compressed())  # two tailed p-value
    return (SCC, P_SCC)

def scipy_pearsonr(Xz, Yz):
    try:
        return scipy.stats.pearsonr(Xz, Yz)
    except:
        return (np.nan, np.nan)

def calc_pearsonr(Xz, Yz):
    mask = np.logical_or(Xz.mask, Yz.mask)
    Xz.mask = mask
    Yz.mask = mask
    (PCC, P_PCC) = scipy_pearsonr(Xz.compressed(), Yz.compressed())  # two tailed p-value
    return (PCC, P_PCC)

def calc_shift_corr(Xz, Yz, D, corfunc=calc_pearsonr):
    d_max = 0
    r_max = 0
    p_max = 1
    for d in range(-D, D+1):
        if d < 0:
            X_seg = Xz[:(len(Xz)+d)]
            Y_seg = Yz[-d:len(Yz)]
        elif d == 0:
            X_seg = Xz
            Y_seg = Yz
        else:
            X_seg = Xz[d:len(Xz)]
            Y_seg = Yz[:len(Yz)-d]
        assert len(X_seg) == len(Y_seg)
        mask = np.logical_or(X_seg.mask, Y_seg.mask)
        X_seg.mask = mask
        Y_seg.mask = mask
        cor = corfunc(X_seg, Y_seg)
        if np.abs(cor[0]) >= np.abs(r_max):
            r_max = cor[0]
            d_max = d
            p_max = cor[1]
    return (r_max, p_max, d_max)

# ... (keep other existing functions)