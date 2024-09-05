import numpy as np
import scipy.stats
from scipy import interpolate
import sys
from .lsalib_core import Q_lam_max, Q_lam_step, my_decimal, Rmax_min, Rmax_max, kcut_min, pipi, pipi_inv
from .lsalib_stats import storeyQvalue, tied_rank, calc_pearsonr, calc_spearmanr, calc_shift_corr, readPvalue, theoPvalue
from .lsalib_normalization import *
from .lsalib_analysis import singleLSA, bootstrapCI, permuPvalue, applyAnalysis
from .lsalib_utils import ma_average, ma_median, sample_wr, fillMissing, safeCmd, float_equal
from . import compcore

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

def storeyQvalue(pvalues, lam=np.arange(0, Q_lam_max, Q_lam_step), method='smoother', robust=False, smooth_df=3):
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

def rpy_spearmanr(Xz, Yz):
    try:
        sr = r('''cor.test''')(Xz, Yz, method='spearman')
        return (sr[3][0], sr[2][0])
    except rpy2.rinterface.RRuntimeError:
        return (np.nan, np.nan)

def rpy_pearsonr(Xz, Yz):
    try:
        sr = r('''cor.test''')(Xz, Yz, method='pearson')
        return (sr[3][0], sr[2][0])
    except rpy2.rinterface.RRuntimeError:
        return (np.nan, np.nan)

def scipy_spearmanr(Xz, Yz):
    try:
        return scipy.stats.spearmanr(Xz, Yz)
    except:
        return (np.nan, np.nan)

def calc_spearmanr(Xz, Yz, sfunc=rpy_spearmanr):
    if not rpy_import:
        sfunc = scipy_spearmanr
    mask = np.logical_or(Xz.mask, Yz.mask)
    Xz.mask = mask
    Yz.mask = mask
    (SCC, P_SCC) = sfunc(Xz.compressed(), Yz.compressed())  # two tailed p-value
    return (SCC, P_SCC)

def scipy_pearsonr(Xz, Yz):
    try:
        return scipy.stats.pearsonr(Xz, Yz)
    except:
        return (np.nan, np.nan)

def calc_pearsonr(Xz, Yz, pfunc=rpy_pearsonr):
    if not rpy_import:
        pfunc = scipy_pearsonr
    mask = np.logical_or(Xz.mask, Yz.mask)
    Xz.mask = mask
    Yz.mask = mask
    (PCC, P_PCC) = pfunc(Xz.compressed(), Yz.compressed())  # two tailed p-value
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

def readPvalue(P_table, R, N, x_sd=1., M=1., alpha=1., beta=1., x_decimal=my_decimal):
    try:
        xi = int(np.around(R*M/(x_sd*np.sqrt(alpha*beta*N))*(10**x_decimal)))
    except OverflowError as ValueError:
        return np.nan
    if xi in P_table:
        return P_table[xi]
    elif xi > max(P_table.keys()):
        return 0.
    else:
        return np.nan

def theoPvalue(Rmax, Dmax=0, precision=.001, x_decimal=my_decimal):
    Rmax = np.max((Rmax, Rmax_min))
    Rmax = np.min((Rmax, Rmax_max))
    print("computing p_table with Rmax=", Rmax, file=sys.stderr)
    P_table = dict()
    for xi in range(0, Rmax*10**(x_decimal)+1):
        if xi == 0:
            P_table[xi] = 1
            continue
        x = xi/float(10**(x_decimal))
        xx = x**2
        pipi_over_xx = pipi/xx
        alpha = precision
        B = 2*Dmax+1
        Kcut = np.max((kcut_min, int(np.ceil(.5 - np.log((alpha/(2**B-1))**(1/B)*xx*(1-np.exp(-pipi_over_xx))/8/2) / pipi_over_xx))))
        A = 1/xx
        Rcdf = 0
        P_two_tail = 1.
        for k in range(1, Kcut+1):
            C = (2*k-1)**2
            Rcdf = Rcdf + (A+pipi_inv/C)*np.exp(-C*pipi_over_xx/2)
            P_current = 1 - (8**B)*(Rcdf**B)
            if P_current < 0:
                break
            else:
                P_two_tail = P_current
        P_table[xi] = P_two_tail
    return P_table

# Add any other statistical functions here