import numpy as np
import scipy as sp
import scipy.linalg
import sys

def calc_tmatrix(bootNum, trend_threshold, timeNum=10000, randomFunc=np.random.normal):
    # return a 3 by 3 transition matrix
    Tm = np.zeros((bootNum, 5))  # each row in order: a,b,c,d,t
    for b in range(0, bootNum):
        Tm[b,] = to_markov(trend_threshold, timeNum, randomFunc)
    return np.average(Tm, axis=0)

def calc_markov_var(P):
    v = 4*((-P[3]/(-1+P[1]+P[2]-2*P[3]))**2)*(1+2*((P[1]-P[2])**2)/(1-((P[1]-P[2])**2)))
    return v

def calc_eigen(P):
    Tmatrix = np.array([P[1], 1-P[1]-P[2], P[2], P[3], 1-2*P[3], P[3], P[2], 1-P[1]-P[2], P[1]]).reshape(3,3)
    w, vl, vr = sp.linalg.eig(Tmatrix, left=True, right=True) 
    w = np.real_if_close(w)
    
    w_one = None
    for i in range(0, len(w)):
        try:
            np.testing.assert_almost_equal(w[i], 1.0)
        except AssertionError:
            continue
        w_one = i
    if w_one is None:
        print("unexpected! reduce float_equal_tolerance", file=sys.stderr)
        quit()
    
    if w_one != 0:  # switch lambda=1 to the first place
        w[0], w[w_one] = w[w_one], w[0]
        vl[:, 0], vl[:, w_one] = vl[:, w_one], vl[:, 0].copy()
        vr[:, 0], vr[:, w_one] = vr[:, w_one], vr[:, 0].copy()
    
    # normalization
    vr[:, 0] = np.array([1.0, 1.0, 1.0])
    vl[:, 0] = vl[:, 0] / np.sum(vl[:, 0])
    for i in range(1, len(w)):
        scale = np.sum(np.inner(vl[:, i], vr[:, i]))
        vl[:, i] = vl[:, i] / scale
    
    # Verification
    try:
        tmp = np.dot(vl.T, vr)
        np.testing.assert_almost_equal(w[0], 1.0)
        np.testing.assert_almost_equal(np.sum(vl[:, 0]), 1.0)
        for i in range(0, len(w)):
            for j in range(0, len(w)):
                if i == j:
                    np.testing.assert_almost_equal(tmp[i, j], 1.0)
                else:
                    np.testing.assert_almost_equal(tmp[i, j], 0.0)
    except AssertionError:
        print("eigen values, right eigen vectors and left eigen vectors not properly found!", file=sys.stderr)
        quit()
    
    return w, vl, vr

def calc_sigma_square(w_sort, vl_sort, vr_sort):
    vl_trans = vl_sort.T
    vr_trans = vr_sort.T
    var_phi = vl_trans[0, :]
    A = (vr_trans[1, 0] - vr_trans[1, 2]) * (vl_trans[1, 0] - vl_trans[1, 2])
    B = (vr_trans[2, 0] - vr_trans[2, 2]) * (vl_trans[2, 0] - vl_trans[2, 2])
    sigma = 4 * var_phi[0]**2 + 2 * var_phi[0]**2 * (
        A**2 * w_sort[1]**2 / (1 - w_sort[1]**2) +
        2 * A * B * w_sort[1] * w_sort[2] / (1 - w_sort[1] * w_sort[2]) +
        B**2 * w_sort[2]**2 / (1 - w_sort[2]**2)
    )
    return sigma

def to_markov(threshold, timeNum, randomFunc):
    t = threshold
    N = timeNum
    Xo = randomFunc(size=N+2)
    Xt = []
    for i in range(0, N+1):
        if Xo[i] == 0 and Xo[i+1] > 0:
            Xt.append(1)
        elif Xo[i] == 0 and Xo[i+1] < 0:
            Xt.append(-1)
        elif Xo[i] == 0 and Xo[i+1] == 0:
            Xt.append(0)
        elif (Xo[i+1] - Xo[i]) / np.abs(Xo[i]) >= t:
            Xt.append(1)
        elif (Xo[i+1] - Xo[i]) / np.abs(Xo[i]) <= -t:
            Xt.append(-1)
        else:
            Xt.append(0)
    
    P = np.zeros(5)
    for j in range(0, N):
        if Xt[j] == 1 and Xt[j+1] == 1:
            P[1] += 1
        elif Xt[j] == 1 and Xt[j+1] == -1:
            P[2] += 1
        elif Xt[j] == 0 and Xt[j+1] == 1:
            P[3] += 1
    
    Xt = np.array(Xt, dtype='int')
    P[0] = np.sum(Xt[Xt == 1])  # di=1
    P[0] = P[0] / (N + 1)
    P[1] = (P[1] / N) / P[0]
    P[2] = (P[2] / N) / P[0]
    P[3] = (P[3] / N) / (1 - 2 * P[0])
    P[4] = t
    return P

def ji_calc_trend(oSeries, lengthSeries, thresh):
    tSeries = np.zeros(lengthSeries-1, dtype='float')
    for i in range(lengthSeries-1):
        if oSeries[i] == 0 and oSeries[i+1] > 0:
            trend = 1
        elif oSeries[i] == 0 and oSeries[i+1] < 0:
            trend = -1
        elif oSeries[i] == 0 and oSeries[i+1] == 0:
            trend = 0
        else:
            trend = (oSeries[i+1]-oSeries[i])/np.abs(oSeries[i])

        if np.isnan(trend):
            tSeries[i] = np.nan
        elif trend >= thresh:
            tSeries[i] = 1
        elif trend <= -thresh:
            tSeries[i] = -1
        else:
            tSeries[i] = 0

    return tSeries