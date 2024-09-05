import numpy as np
import scipy as sp
from .lsalib_stats import storeyQvalue, tied_rank, calc_pearsonr, calc_spearmanr, calc_shift_corr, readPvalue
from .lsalib_utils import ma_average, sample_wr, fillMissing
from .lsalib_normalization import *
from .lsalib_trend import ji_calc_trend
from .lsalib_core import P_table
from . import compcore

def singleLSA(series1, series2, delayLimit, fTransform, zNormalize, trendThresh=None, keepTrace=True):
    timespots = series1.shape[1]
    if trendThresh is not None:
        xSeries = ji_calc_trend(zNormalize(fTransform(series1)), timespots, trendThresh)
        ySeries = ji_calc_trend(zNormalize(fTransform(series2)), timespots, trendThresh)
    else:
        xSeries = zNormalize(fTransform(series1))
        ySeries = zNormalize(fTransform(series2))

    lsad = compcore.LSA_Data(delayLimit, xSeries, ySeries)
    lsar = compcore.DP_lsa(lsad, keepTrace)
    del lsad
    return lsar

def bootstrapCI(series1, series2, Smax, delayLimit, bootCI, bootNum, fTransform, zNormalize, trendThresh=None, debug=0):
    if series1.shape[0] == 1:
        return (Smax, Smax, Smax)

    assert series1.shape[1] == series2.shape[1]
    timespots = series1.shape[1]
    lengthSeries = series1.shape[1]
    if trendThresh is not None:
        lengthSeries = timespots - 1

    lsad = compcore.LSA_Data()
    BS_set = np.zeros(bootNum, dtype='float')
    for i in range(bootNum):
        if trendThresh is None:
            Xb = zNormalize(fTransform(np.ma.array([sample_wr(series1[:,j], series1.shape[0]) for j in range(series1.shape[1])]).T))
            Yb = zNormalize(fTransform(np.ma.array([sample_wr(series2[:,j], series2.shape[0]) for j in range(series2.shape[1])]).T))
        else:
            Xb = ji_calc_trend(zNormalize(fTransform(np.ma.array([sample_wr(series1[:,j], series1.shape[0]) for j in range(series1.shape[1])]).T)), timespots, trendThresh)
            Yb = ji_calc_trend(zNormalize(fTransform(np.ma.array([sample_wr(series2[:,j], series2.shape[0]) for j in range(series2.shape[1])]).T)), timespots, trendThresh)
        
        lsad.assign(delayLimit, Xb, Yb)
        BS_set[i] = compcore.DP_lsa(lsad, False).score

    BS_set.sort()
    BS_mean = np.mean(BS_set)
    a1 = (1-bootCI)/2.0
    a2 = bootCI+(1-bootCI)/2.0

    if debug in [1, 3]:
        Smax = 2*Smax - BS_mean
    if debug in [2, 3]:
        z0 = sp.stats.norm.ppf(np.sum(BS_set <= Smax)/float(bootNum))
        a1 = sp.stats.norm.cdf(2*z0+sp.stats.norm.ppf(a1))
        a2 = sp.stats.norm.cdf(2*z0+sp.stats.norm.ppf(a2))

    return (BS_mean, BS_set[int(np.floor(bootNum*a1))-1], BS_set[int(np.ceil(bootNum*a2))-1])

def permuPvalue(series1, series2, delayLimit, precisionP, Smax, fTransform, zNormalize, trendThresh=None):
    lengthSeries = series1.shape[1]
    timespots = series1.shape[1]
    if trendThresh is not None:
        lengthSeries = timespots - 1

    lsad = compcore.LSA_Data()
    PP_set = np.zeros(precisionP, dtype='float')
    if trendThresh is None:
        Xz = zNormalize(fTransform(series1))
    else:
        Xz = ji_calc_trend(zNormalize(fTransform(series1)), timespots, trendThresh)
    Y = np.ma.array(series2)

    for i in range(precisionP):
        np.random.shuffle(Y.T)
        if trendThresh is None:
            Yz = zNormalize(fTransform(Y))
        else:
            Yz = ji_calc_trend(zNormalize(fTransform(Y)), timespots, trendThresh)
        lsad.assign(delayLimit, Xz, Yz)
        PP_set[i] = compcore.DP_lsa(lsad, False).score

    if Smax >= 0:
        P_two_tail = np.sum(np.abs(PP_set) >= Smax)/float(precisionP)
    else:
        P_two_tail = np.sum(-np.abs(PP_set) <= Smax)/float(precisionP)
    return P_two_tail

def applyAnalysis(firstData, secondData, onDiag=True, delayLimit=3, minOccur=.5, bootCI=.95, bootNum=0, pvalueMethod='perm', precisionP=1000,
                  fTransform=None, zNormalize=None, approxVar=1, resultFile=None, trendThresh=None,
                  firstFactorLabels=None, secondFactorLabels=None, qvalueMethod='scipy', progressive=0):
    firstFactorNum, secondFactorNum = firstData.shape[0], secondData.shape[0]
    spotNum = firstData.shape[2]
    
    lsaTable = [None] * (firstFactorNum * secondFactorNum)
    pvalues = np.zeros(firstFactorNum * secondFactorNum)
    pccpvalues = np.zeros(firstFactorNum * secondFactorNum)
    spccpvalues = np.zeros(firstFactorNum * secondFactorNum)
    sccpvalues = np.zeros(firstFactorNum * secondFactorNum)
    ssccpvalues = np.zeros(firstFactorNum * secondFactorNum)
    
    ti = 0
    for i in range(firstFactorNum):
        for j in range(secondFactorNum):
            if onDiag and j <= i:
                continue
            
            Xz = firstData[i]
            Yz = secondData[j]
            
            if np.all(np.isnan(Xz)) or np.all(np.isnan(Yz)):
                lsaTable[ti] = [i, j, 0, 0, 0, -1, -1, 0, 0, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                pvalues[ti] = np.nan
                pccpvalues[ti] = np.nan
                spccpvalues[ti] = np.nan
                sccpvalues[ti] = np.nan
                ssccpvalues[ti] = np.nan
            else:
                LSA_result = singleLSA(Xz, Yz, delayLimit, fTransform, zNormalize, trendThresh, True)
                
                Smax = LSA_result.score
                Al = len(LSA_result.trace)
                PCC, P_PCC = calc_pearsonr(ma_average(Xz, axis=0), ma_average(Yz, axis=0))
                SCC, P_SCC = calc_spearmanr(ma_average(Xz, axis=0), ma_average(Yz, axis=0))
                
                try:
                    SPCC, P_SPCC, D_SPCC = calc_shift_corr(ma_average(Xz, axis=0), ma_average(Yz, axis=0), delayLimit, calc_pearsonr)
                    SSCC, P_SSCC, D_SSCC = calc_shift_corr(ma_average(Xz, axis=0), ma_average(Yz, axis=0), delayLimit, calc_spearmanr)
                except FloatingPointError:
                    SPCC, P_SPCC, D_SPCC = np.nan, np.nan, np.nan
                    SSCC, P_SSCC, D_SSCC = np.nan, np.nan, np.nan
                
                pccpvalues[ti] = P_PCC
                spccpvalues[ti] = P_SPCC
                sccpvalues[ti] = P_SCC
                ssccpvalues[ti] = P_SSCC
                
                if Al == 0:
                    pvalues[ti] = np.nan
                    lsaTable[ti] = [i, j, 0, 0, 0, -1, -1, 0, 0, 1, PCC, P_PCC, SPCC, P_SPCC, D_SPCC, SCC, P_SCC, SSCC, P_SSCC, D_SSCC]
                else:
                    Xs, Ys, Al = LSA_result.trace[Al-1][0], LSA_result.trace[Al-1][1], len(LSA_result.trace)
                    
                    lsaP = -1
                    if pvalueMethod in ['theo', 'mix']:
                        lsaP = readPvalue(P_table, R=np.abs(Smax)*spotNum, N=spotNum, x_sd=1, M=1, alpha=1., beta=1., x_decimal=2)
                    
                    if (pvalueMethod in ['mix'] and lsaP <= 0.05) or (pvalueMethod in ['perm']):
                        Xp = np.ma.array(Xz, copy=True)
                        Yp = np.ma.array(Yz, copy=True)
                        lsaP = permuPvalue(Xp, Yp, delayLimit, precisionP, np.abs(Smax), fTransform, zNormalize, trendThresh)
                    
                    pvalues[ti] = lsaP
                    
                    if bootNum > 0:
                        Xb = np.ma.array(Xz, copy=True)
                        Yb = np.ma.array(Yz, copy=True)
                        Smax, Sl, Su = bootstrapCI(Xb, Yb, Smax, delayLimit, bootCI, bootNum, fTransform, zNormalize, trendThresh)
                    else:
                        Sl, Su = Smax, Smax
                    
                    lsaTable[ti] = [i, j, Smax, Sl, Su, Xs, Ys, Al, Xs-Ys, lsaP, PCC, P_PCC, SPCC, P_SPCC, D_SPCC, SCC, P_SCC, SSCC, P_SSCC, D_SSCC]
            
            ti += 1
    
    qvalues = storeyQvalue(pvalues)
    pccqvalues = storeyQvalue(pccpvalues)
    spccqvalues = storeyQvalue(spccpvalues)
    sccqvalues = storeyQvalue(sccpvalues)
    ssccqvalues = storeyQvalue(ssccpvalues)
    
    for k in range(len(qvalues)):
        lsaTable[k].extend([qvalues[k], pccqvalues[k], spccqvalues[k], sccqvalues[k], ssccqvalues[k]])
    
    if resultFile:
        for row in lsaTable:
            if row:
                print("\t".join(['%s']*len(row)) % tuple(row), file=resultFile)
    
    return lsaTable