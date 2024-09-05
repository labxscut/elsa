import numpy as np
import scipy as sp
from .lsalib_stats import storeyQvalue, tied_rank, calc_pearsonr, calc_spearmanr, calc_shift_corr
from .lsalib_stats import ma_average, sample_wr, fillMissing
from .lsalib_normalization import *
from .lsalib_trend import ji_calc_trend
from .lsalib_core import P_table
from .lsalib_theory import readPvalue, theoPvalue
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
    timespots = series1.shape[1]
    if trendThresh is not None:
        xSeries = ji_calc_trend(zNormalize(fTransform(series1)), timespots, trendThresh)
        ySeries = ji_calc_trend(zNormalize(fTransform(series2)), timespots, trendThresh)
    else:
        xSeries = zNormalize(fTransform(series1))
        ySeries = zNormalize(fTransform(series2))

    bootResults = []
    for _ in range(bootNum):
        bootIndex = sample_wr(range(xSeries.shape[0]), xSeries.shape[0])
        xBoot = xSeries[bootIndex, :]
        yBoot = ySeries[bootIndex, :]
        
        lsad = compcore.LSA_Data(delayLimit, xBoot, yBoot)
        lsar = compcore.DP_lsa(lsad, False)
        bootResults.append(lsar.score)

    bootResults = np.array(bootResults)
    ciLow = np.percentile(bootResults, (1 - bootCI) * 100 / 2)
    ciHigh = np.percentile(bootResults, (1 + bootCI) * 100 / 2)

    return ciLow, ciHigh

def permuPvalue(series1, series2, delayLimit, precisionP, Smax, fTransform, zNormalize, trendThresh=None):
    timespots = series1.shape[1]
    if trendThresh is not None:
        xSeries = ji_calc_trend(zNormalize(fTransform(series1)), timespots, trendThresh)
        ySeries = ji_calc_trend(zNormalize(fTransform(series2)), timespots, trendThresh)
    else:
        xSeries = zNormalize(fTransform(series1))
        ySeries = zNormalize(fTransform(series2))

    lsad = compcore.LSA_Data(delayLimit, xSeries, ySeries)
    lsar = compcore.DP_lsa(lsad, False)
    
    permResults = []
    for _ in range(precisionP):
        permIndex = np.random.permutation(xSeries.shape[0])
        xPerm = xSeries[permIndex, :]
        yPerm = ySeries
        
        lsad_perm = compcore.LSA_Data(delayLimit, xPerm, yPerm)
        lsar_perm = compcore.DP_lsa(lsad_perm, False)
        permResults.append(lsar_perm.score)

    permResults = np.array(permResults)
    pvalue = np.sum(permResults >= lsar.score) / precisionP

    return pvalue

def applyAnalysis(firstData, secondData, onDiag=True, delayLimit=3, minOccur=.5, bootCI=.95, bootNum=0, pvalueMethod='perm', precisionP=1000,
                  fTransform=None, zNormalize=None, approxVar=1, resultFile=None, trendThresh=None,
                  firstFactorLabels=None, secondFactorLabels=None, qvalueMethod='scipy', progressive=0):
    # Initialize variables
    resultList = []
    totalPairs = len(firstData) * (len(secondData) if not onDiag else 1)
    pairsDone = 0

    # Iterate through data pairs
    for i, first in enumerate(firstData):
        for j, second in enumerate(secondData):
            if onDiag and i != j:
                continue

            # Calculate LSA score
            lsad = compcore.LSA_Data(delayLimit, first, second)
            lsar = compcore.DP_lsa(lsad, False)
            
            # Calculate confidence interval if bootNum > 0
            if bootNum > 0:
                ciLow, ciHigh = bootstrapCI(first, second, delayLimit, bootCI, bootNum, fTransform, zNormalize, trendThresh)
            else:
                ciLow, ciHigh = None, None

            # Calculate p-value
            if pvalueMethod == 'perm':
                pvalue = permuPvalue(first, second, delayLimit, precisionP, lsar.score, fTransform, zNormalize, trendThresh)
            elif pvalueMethod == 'theo':
                pvalue = readPvalue(P_table, lsar.score, first.shape[1], x_sd=1., M=1., alpha=1., beta=1.)
            else:
                raise ValueError("Invalid pvalueMethod. Choose 'perm' or 'theo'.")

            # Append results
            resultList.append({
                'x': i, 'y': j,
                'score': lsar.score, 'delay': lsar.delay,
                'pvalue': pvalue, 'ciLow': ciLow, 'ciHigh': ciHigh
            })

            # Update progress
            pairsDone += 1
            if progressive and pairsDone % progressive == 0:
                print(f"Progress: {pairsDone}/{totalPairs} pairs processed", file=sys.stderr)

    # Convert results to DataFrame
    results = pd.DataFrame(resultList)

    # Calculate q-values
    if qvalueMethod == 'scipy':
        _, results['qvalue'], _, _ = multipletests(results['pvalue'], method='fdr_bh')
    elif qvalueMethod == 'storey':
        results['qvalue'] = storeyQvalue(results['pvalue'])
    else:
        raise ValueError("Invalid qvalueMethod. Choose 'scipy' or 'storey'.")

    # Save results if resultFile is provided
    if resultFile:
        results.to_csv(resultFile, index=False)

    return results

# Add other analysis functions here