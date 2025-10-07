#lalib.py -- Library of Liquid Association Analysis(LAA) Package
#LICENSE: BSD

#Copyright (c) 2008 Li Charles Xia
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions
#are met:
#1. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#3. The name of the author may not be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
#IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
#IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
#NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""llalib.py -- Library of Liquid Association Analysis(LA) Package

  NOTE: numpy and scipy is required to use thie module
  NOTE: accepts input sequence table as delimited text file 
        with first row is the "Date" and other factor labels and
        first column is the time spot labels
"""

import csv
import sys
import os
import random
import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.stats

try:
    # When running as installed package
    from lsa import lsalib
    from lsa import compcore
except ImportError:
    # When running in debug mode
    try:
        from . import lsalib
        from . import compcore
    except ImportError:
        import lsalib
        import compcore

#global variable, stores calculated p-values.
disp_decimal = 8
kcut_min = 100
Rmax_min = 10
my_decimal = 2    # preset x step size for P_table
pipi = np.pi**2  # pi^2
pipi_inv = 1/pipi
Q_lam_step = 0.05
Q_lam_max = 0.95

### LA with scoutVars functions ###

def applyLA(inputData, scoutVars, factorLabels, bootCI=.95, bootNum=1000, minOccur=.50, 
           pvalueMethod=1000, fTransform=lsalib.simpleAverage, 
           zNormalize=lsalib.noZeroNormalize, resultFile=None):

    col_labels = ['X','Y','Z','LA','lowCI','upCI','P','Q','Xi','Yi','Zi']
    print("\t".join(col_labels), file=resultFile)

    inputFactorNum = inputData.shape[0]
    inputRepNum = inputData.shape[1]
    inputSpotNum = inputData.shape[2]
    scoutNum = len(scoutVars)
    cp = np.array([False]*inputFactorNum**3, dtype='bool') #consider bitvector
    cp.shape = (inputFactorNum, inputFactorNum, inputFactorNum)
    laTable = [None]*inputFactorNum*scoutNum
    pvalues = np.zeros(inputFactorNum*scoutNum, dtype='float')
    timespots = inputSpotNum #same length already assumed
    replicates = inputRepNum
    ti = 0
    
    for i in range(0, scoutNum):
        Xi = scoutVars[i][0] - 1
        Yi = scoutVars[i][1] - 1
        Xo = np.ma.masked_invalid(inputData[Xi], copy=True)
        Yo = np.ma.masked_invalid(inputData[Yi], copy=True)
        
        for j in range(0, inputFactorNum):
            Zi = j
            if Xi == Yi or Xi == Zi or Zi == Yi:
                continue   #ignore invalid entries
            if cp[Xi,Yi,Zi] or cp[Xi,Zi,Yi] or cp[Zi,Xi,Yi] or cp[Zi,Yi,Xi] or cp[Yi,Xi,Zi] or cp[Yi,Zi,Xi]:
                continue   #ignore redundant entries
                
            cp[Xi,Yi,Zi] = True
            Zo = np.ma.masked_invalid(inputData[Zi], copy=True)
            
            Xo_minOccur = np.sum(np.logical_or(np.isnan(lsalib.ma_average(Xo)), 
                                lsalib.ma_average(Xo)==0))/float(timespots) < minOccur
            Yo_minOccur = np.sum(np.logical_or(np.isnan(lsalib.ma_average(Yo)), 
                                lsalib.ma_average(Yo)==0))/float(timespots) < minOccur
            Zo_minOccur = np.sum(np.logical_or(np.isnan(lsalib.ma_average(Zo)), 
                                lsalib.ma_average(Zo)==0))/float(timespots) < minOccur
                                
            if Xo_minOccur or Yo_minOccur or Zo_minOccur:
                continue
            if np.all(Xo.mask) or np.all(Yo.mask) or np.all(Zo.mask):
                continue
                
            LA_score = singleLA(Xo, Yo, Zo, fTransform, zNormalize)

            if pvalueMethod >= 0 or pvalueMethod < 0:
                Xp = np.ma.array(Xo, copy=True)
                Yp = np.ma.array(Yo, copy=True)
                Zp = np.ma.array(Zo, copy=True)
                laP = LApermuPvalue(Xp, Yp, Zp, abs(pvalueMethod), np.abs(LA_score), 
                                  fTransform, zNormalize)
            else:
                print("This should not be reached now", file=sys.stderr)
                
            pvalues[ti] = laP
            
            if bootNum > 0:
                Xb = np.ma.array(Xo, copy=True)
                Yb = np.ma.array(Yo, copy=True)
                Zb = np.ma.array(Zo, copy=True)
                (LA_score, Sl, Su) = LAbootstrapCI(Xb, Yb, Zb, LA_score, bootCI, bootNum, 
                                                 fTransform, zNormalize)
            else:
                (LA_score, Sl, Su) = (LA_score, LA_score, LA_score)

            laTable[ti] = [Xi, Yi, Zi, LA_score, Sl, Su, laP]
            ti += 1

    pvalues = pvalues[:ti]
    laTable = laTable[:ti]
    qvalues = lsalib.storeyQvalue(pvalues)
    
    for k in range(0, len(qvalues)):
        laTable[k] = laTable[k] + [qvalues[k], laTable[k][0]+1, laTable[k][1]+1, laTable[k][2]+1]

    for row in laTable:
        print("\t".join(['%s']*len(col_labels)) % 
              tuple([factorLabels[row[0]], factorLabels[row[1]], factorLabels[row[2]]] + 
                    [f"{v:.{disp_decimal}f}" if isinstance(v, float) else v for v in row[3:]]), 
              file=resultFile)

def calc_LA(series1, series2, series3):     # the function to calculate LA score
    n1 = len(series1)
    n2 = len(series2)
    n3 = len(series3)
    assert n1==n2 and n2 == n3
    return np.sum(series1*series2*series3)/n1

def singleLA(series1, series2, series3, fTransform, zNormalize): # calculate normalized LA score
    return calc_LA(zNormalize(fTransform(series1)),zNormalize(fTransform(series2)),zNormalize(fTransform(series3)))

def LAbootstrapCI(series1, series2, series3, LA_score, bootCI, bootNum, fTransform, zNormalize, debug=0):
    ### no feasible, skipping bootstraping
    if series1.shape[0] == 1:
        return (LA_score, LA_score, LA_score)

    BS_set = np.zeros(bootNum, dtype='float')
    for i in range(0, bootNum):
        Xb = np.ma.array([ lsalib.sample_wr(series1[:,j], series1.shape[0]) for j in range(0,series1.shape[1]) ]).T
        Yb = np.ma.array([ lsalib.sample_wr(series2[:,j], series2.shape[0]) for j in range(0,series2.shape[1]) ]).T
        Zb = np.ma.array([ lsalib.sample_wr(series3[:,j], series3.shape[0]) for j in range(0,series3.shape[1]) ]).T
        BS_set[i] = calc_LA(Xb, Yb, Zb)
    BS_set.sort()                                 #from smallest to largest
    BS_mean = np.mean(BS_set)
    return ( BS_mean, BS_set[np.floor(bootNum*a1)-1], BS_set[np.ceil(bootNum*a2)-1] )

def LApermuPvalue(series1, series2, series3, pvalueMethod, LA_score, fTransform, zNormalize):
    PP_set = np.zeros(pvalueMethod, dtype='float')
    X = zNormalize(fTransform(series1))
    Y = zNormalize(fTransform(series2))
    Z = np.ma.array(series3)                                               #use = only assigns reference, must use a constructor
    for i in range(0, pvalueMethod):
        np.random.shuffle(Z.T)
        PP_set[i] = calc_LA(X, Y, zNormalize(fTransform(Z)))
    if LA_score >= 0:
        P_two_tail = np.sum(np.abs(PP_set) >= LA_score)/float(pvalueMethod)
    else:
        P_two_tail = np.sum(-np.abs(PP_set) <= LA_score)/float(pvalueMethod)
    return P_two_tail


# ab initio local liquid association analysis functions ###  

def transform_series(series, fTransform, zNormalize):
    """Transform a single time series while preserving masked array structure."""
    # Ensure input is masked array
    if not isinstance(series, np.ma.MaskedArray):
        series = np.ma.array(series)
    
    # Apply transform while preserving mask
    transformed = fTransform(series)
    if not isinstance(transformed, np.ma.MaskedArray):
        transformed = np.ma.array(transformed, mask=series.mask)
    
    # Apply normalization while preserving mask
    normalized = zNormalize(transformed)
    if not isinstance(normalized, np.ma.MaskedArray):
        normalized = np.ma.array(normalized, mask=series.mask)
    
    return normalized

def applyLLAnalysis(cleanData, factorLabels, delayLimit=3, bootCI=.95, bootNum=1000, minOccur=.50,
                   pvalueMethod="perm", precision=1000, fillMethod='linear', normMethod='pnz',
                   fTransform=lsalib.simpleAverage, zNormalize=lsalib.noZeroNormalize, 
                   resultFile=None, qvalue_func=lsalib.storeyQvalue):
    """Apply Local Liquid Association analysis to input data."""
    col_labels = ['X','Y','Z','LA','lowCI','upCI','P','Q','Xi','Yi','Zi','Delay']
    print("\t".join(col_labels), file=resultFile)
    
    inputFactorNum = cleanData.shape[0]
    inputRepNum = cleanData.shape[1]
    inputSpotNum = cleanData.shape[2]
    
    laTable = []
    pvalues = []
    
    # Track processed triplets to avoid redundancy
    processed = set()
    
    for Xi in range(inputFactorNum):
        for Yi in range(Xi + 1, inputFactorNum):
            for Zi in range(Yi + 1, inputFactorNum):
                triplet = tuple(sorted([Xi, Yi, Zi]))
                if triplet in processed:
                    continue
                processed.add(triplet)
                
                try:
                    # Transform data while preserving masked array structure
                    X = transform_series(cleanData[Xi], fTransform, zNormalize)
                    Y = transform_series(cleanData[Yi], fTransform, zNormalize)
                    Z = transform_series(cleanData[Zi], fTransform, zNormalize)
                    
                    # Check minimum occurrence criteria
                    if not all(np.sum(~np.ma.getmask(s))/float(inputSpotNum) >= minOccur 
                             for s in [X, Y, Z]):
                        continue
                    
                    # Calculate LA score and statistics
                    lla_data = compcore.LLA_Data(delayLimit, X, Y, Z)
                    lla_result = compcore.DP_lla(lla_data)
                    
                    # Calculate p-value
                    pvalue = (LLApermuPvalue(X, Y, Z, delayLimit, precision, lla_result.score) 
                            if pvalueMethod == "perm" else pvalueMethod)
                    pvalues.append(pvalue)
                    
                    # Calculate bootstrap CI if requested
                    if bootNum > 0:
                        la_score, lowCI, upCI = LLAbootstrapCI(X, Y, Z, lla_result.score, 
                                                           delayLimit, bootCI, bootNum)
                    else:
                        la_score = lowCI = upCI = lla_result.score
                    
                    laTable.append([Xi, Yi, Zi, la_score, lowCI, upCI, pvalue])
                    
                except Exception as e:
                    print("Error during analysis:", file=sys.stderr)
                    traceback.print_exc()
                    continue
    
    # Calculate q-values and write results
    if pvalues:
        qvalues = qvalue_func(np.array(pvalues))
        for k, row in enumerate(laTable):
            laTable[k] = row[:7] + [qvalues[k]] + [x+1 for x in row[:3]] + [0]
            
        for row in laTable:
            print("\t".join(['%s']*len(col_labels)) % 
                  tuple([factorLabels[row[0]], factorLabels[row[1]], factorLabels[row[2]]] + 
                        [f"{v:.8f}" if isinstance(v, float) else v for v in row[3:]]), 
                  file=resultFile)
    else:
        print("No valid triplets found for analysis", file=sys.stderr)

def LLApermuPvalue(X, Y, Z, delayLimit, precisionP, LA_score):
    """Compute permutation-based p-value for LLA score.
    
    Args:
        X, Y, Z: Masked arrays containing time series data
        delayLimit: Maximum time delay to consider
        precisionP: Number of permutations
        LA_score: Observed LA score
        
    Returns:
        float: Two-tailed p-value
    """
    PP_set = np.zeros(precisionP, dtype='float')
    
    # Ensure inputs are masked arrays
    X = np.ma.array(X) if not isinstance(X, np.ma.MaskedArray) else X
    Y = np.ma.array(Y) if not isinstance(Y, np.ma.MaskedArray) else Y
    Z = np.ma.array(Z) if not isinstance(Z, np.ma.MaskedArray) else Z

    for i in range(precisionP):
        # Shuffle Z while preserving mask
        Zp = np.ma.array(Z)
        idx = np.random.permutation(len(Z))
        Zp = Zp[idx]
        
        lla_data = compcore.LLA_Data(delayLimit, X, Y, Zp)
        PP_set[i] = compcore.DP_lla(lla_data).score

    if LA_score >= 0:
        P_two_tail = np.sum(np.abs(PP_set) >= LA_score) / float(precisionP)
    else:
        P_two_tail = np.sum(-np.abs(PP_set) <= LA_score) / float(precisionP)
    return P_two_tail

def LLAbootstrapCI(X, Y, Z, LA_score, delayLimit, bootCI, bootNum):
    """Compute bootstrap confidence intervals for LLA score.
    
    Args:
        X, Y, Z: Masked arrays containing time series data
        LA_score: Original LA score
        delayLimit: Maximum time delay
        bootCI: Confidence interval (e.g., 0.95)
        bootNum: Number of bootstrap iterations
        
    Returns:
        tuple: (mean LA score, lower CI, upper CI)
    """
    if len(X) <= 1:  # Not enough data points
        return LA_score, LA_score, LA_score
        
    BS_set = np.zeros(bootNum, dtype='float')
    
    # Ensure inputs are masked arrays
    X = np.ma.array(X) if not isinstance(X, np.ma.MaskedArray) else X
    Y = np.ma.array(Y) if not isinstance(Y, np.ma.MaskedArray) else Y
    Z = np.ma.array(Z) if not isinstance(Z, np.ma.MaskedArray) else Z
    
    for i in range(bootNum):
        # Bootstrap sampling with replacement
        idx = np.random.choice(len(X), size=len(X), replace=True)
        Xb = X[idx]
        Yb = Y[idx]
        Zb = Z[idx]
        
        lla_data = compcore.LLA_Data(delayLimit, Xb, Yb, Zb)
        BS_set[i] = compcore.DP_lla(lla_data).score
        
    BS_set.sort()
    BS_mean = np.mean(BS_set)
    
    # Calculate CI indices
    alpha = 1 - bootCI
    lower_idx = int(np.floor(bootNum * (alpha/2)))
    upper_idx = int(np.ceil(bootNum * (1 - alpha/2))) - 1
    
    return BS_mean, BS_set[lower_idx], BS_set[upper_idx]