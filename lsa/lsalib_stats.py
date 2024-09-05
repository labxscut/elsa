import numpy as np
import scipy.stats
from scipy import interpolate

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

# Add other statistical functions here