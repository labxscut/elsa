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
import json
import random
import traceback
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
        from lsa import lsalib
        from lsa import compcore

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
    # Calculate CI indices based on bootCI
    alpha = 1 - bootCI
    lower_idx = int(np.floor(bootNum * (alpha/2)))
    upper_idx = int(np.ceil(bootNum * (1 - alpha/2))) - 1
    return ( BS_mean, BS_set[lower_idx], BS_set[upper_idx] )

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

def applyLLAnalysis(cleanData, factorLabels, delayLimit, bootCI=.95, bootNum=1000, minOccur=.50,
                   pvalueMethod="perm", precision=1000, fillMethod='linear', normMethod='pnz',
                   fTransform=lsalib.simpleAverage, zNormalize=lsalib.noZeroNormalize, 
                   resultFile=None, qvalue_func=lsalib.storeyQvalue, keep_trace=False,
                   checkpoint_file=None, resume=False, flush_every=100, checkpoint_every=1000):
    """Apply Local Liquid Association analysis to input data."""
    # Define column format specifications
    col_formats = {
        'X': ('%-6s', '%-6s'),      # (header_fmt, data_fmt) - Name columns left-aligned
        'Y': ('%-6s', '%-6s'),
        'Z': ('%-6s', '%-6s'),
        'LLA': ('%-8s', '%-8.4f'),     # Numeric columns with 4 decimal places
        'Delay': ('%-6s', '%-6d'),     # Y-Z delay (X-Y always synchronous)
        'Start_X': ('%-8s', '%-8d'),
        'Start_Y': ('%-8s', '%-8d'),
        'Start_Z': ('%-8s', '%-8d'),
        'End_X': ('%-6s', '%-6d'),
        'End_Y': ('%-6s', '%-6d'),
        'End_Z': ('%-6s', '%-6d'),
        'lowCI': ('%-8s', '%-8.4f'),
        'upCI': ('%-8s', '%-8.4f'),
        'P': ('%-8s', '%-8.4f')
    }
    
    # Column order
    columns = ['X', 'Y', 'Z', 'LLA', 
              'Start_X', 'Start_Y', 'Start_Z', 'End_X', 'End_Y', 'End_Z',
              'P', 'lowCI', 'upCI', 'Delay']
    
    # Skip header when appending to an existing non-empty result file
    file_has_content = False
    if hasattr(resultFile, 'tell') and hasattr(resultFile, 'seek'):
        cur_pos = resultFile.tell()
        resultFile.seek(0, os.SEEK_END)
        file_has_content = resultFile.tell() > 0
        resultFile.seek(cur_pos, os.SEEK_SET)

    if not (resume and file_has_content):
        header_parts = [col_formats[col][0] % col for col in columns]
        print(' '.join(header_parts), file=resultFile)
        resultFile.flush()
    
    inputFactorNum = cleanData.shape[0]
    inputRepNum = cleanData.shape[1]
    inputSpotNum = cleanData.shape[2]
    
    # Restore checkpoint if requested
    resume_from = (0, 1, 2)
    completed_rows = 0
    if resume and checkpoint_file and os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, 'r') as ckf:
                state = json.load(ckf)
            if state.get('completed', False):
                print("Checkpoint indicates completed analysis; skipping computation.", file=sys.stderr)
                return
            last_triplet = state.get('last_triplet')
            if isinstance(last_triplet, list) and len(last_triplet) == 3:
                lx, ly, lz = [int(v) for v in last_triplet]
                if lz + 1 < inputFactorNum:
                    resume_from = (lx, ly, lz + 1)
                elif ly + 1 < inputFactorNum - 1:
                    resume_from = (lx, ly + 1, ly + 2)
                elif lx + 1 < inputFactorNum - 2:
                    resume_from = (lx + 1, lx + 2, lx + 3)
                else:
                    print("Checkpoint points to end of search space; marking complete.", file=sys.stderr)
                    return
            completed_rows = int(state.get('rows_written', 0))
            print(f"Resuming from triplet index {resume_from}, rows_written={completed_rows}", file=sys.stderr)
        except Exception:
            print("Warning: failed to load checkpoint; restarting from beginning.", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            resume_from = (0, 1, 2)

    def _before_resume_point(xi, yi, zi, start):
        return (xi < start[0] or
                (xi == start[0] and yi < start[1]) or
                (xi == start[0] and yi == start[1] and zi < start[2]))

    def _write_checkpoint(last_triplet, rows_written, completed=False):
        if not checkpoint_file:
            return
        state = {
            'last_triplet': list(last_triplet) if last_triplet is not None else None,
            'rows_written': int(rows_written),
            'completed': bool(completed)
        }
        tmp_path = checkpoint_file + '.tmp'
        with open(tmp_path, 'w') as ckf:
            json.dump(state, ckf)
            ckf.flush()
            os.fsync(ckf.fileno())
        os.replace(tmp_path, checkpoint_file)

    # Precompute transformed series and occurrence validity once per factor
    transformed = [None] * inputFactorNum
    valid_occurrence = [False] * inputFactorNum
    for i in range(inputFactorNum):
        series_i = transform_series(cleanData[i], fTransform, zNormalize)
        transformed[i] = series_i
        valid_occurrence[i] = (np.sum(~np.ma.getmask(series_i)) / float(inputSpotNum)) >= minOccur

    written_since_flush = 0
    last_triplet_written = None
    
    for Xi in range(inputFactorNum):
        for Yi in range(Xi + 1, inputFactorNum):
            for Zi in range(Yi + 1, inputFactorNum):
                if resume and _before_resume_point(Xi, Yi, Zi, resume_from):
                    continue
                
                try:
                    X = transformed[Xi]
                    Y = transformed[Yi]
                    Z = transformed[Zi]
                    
                    # Check minimum occurrence criteria
                    if not (valid_occurrence[Xi] and valid_occurrence[Yi] and valid_occurrence[Zi]):
                        continue
                    
                    # Calculate LA score and statistics
                    lla_data = compcore.LLA_Data(delayLimit, X, Y, Z)
                    lla_result = compcore.DP_lla(lla_data, keep_trace=keep_trace) # keep_trace
                    
                    # 默认值（当 keep_trace=False 时）
                    start_triplet = [-1, -1, -1]
                    end_triplet = [-1, -1, -1]

                    if keep_trace and lla_result.trace:
                        end_triplet = lla_result.trace[0]
                        start_triplet = lla_result.trace[-1]
                        # C++ trace 索引为 1-based，可直接用于输出
                    
                    # Calculate p-value
                    if pvalueMethod == "perm":
                        pvalue = LLApermuPvalue(X, Y, Z, delayLimit, precision, lla_result.score)
                    elif pvalueMethod == "theo":
                        pvalue = LLAtheoPvalue(lla_result.score, inputSpotNum, delayLimit, precision)
                    elif pvalueMethod == "mix":
                        theo_p = LLAtheoPvalue(lla_result.score, inputSpotNum, delayLimit, precision)
                        if theo_p <= 0.05:
                            pvalue = LLApermuPvalue(X, Y, Z, delayLimit, precision, lla_result.score)
                        else:
                            pvalue = theo_p
                    else:
                        pvalue = float(pvalueMethod)
                    
                    # Calculate bootstrap CI if requested
                    if bootNum > 0:
                        la_score, lowCI, upCI = LLAbootstrapCI(X, Y, Z, lla_result.score, 
                                                           delayLimit, bootCI, bootNum)
                    else:
                        la_score = lowCI = upCI = lla_result.score
                    
                    # Calculate Delay (Y-Z only, since X-Y are always synchronous)
                    # Note: X-Y are always synchronous (i==j), so we only track Y-Z delay
                    if keep_trace and (start_triplet[0] != -1) and (end_triplet[0] != -1):
                        delay = end_triplet[2] - end_triplet[1]  # Y-Z delay (can be non-zero)
                    else:
                        delay = 0

                    data_values = [
                        col_formats['X'][1] % factorLabels[Xi],
                        col_formats['Y'][1] % factorLabels[Yi],
                        col_formats['Z'][1] % factorLabels[Zi],
                        col_formats['LLA'][1] % la_score,
                        col_formats['Start_X'][1] % start_triplet[0],
                        col_formats['Start_Y'][1] % start_triplet[1],
                        col_formats['Start_Z'][1] % start_triplet[2],
                        col_formats['End_X'][1] % end_triplet[0],
                        col_formats['End_Y'][1] % end_triplet[1],
                        col_formats['End_Z'][1] % end_triplet[2],
                        col_formats['P'][1] % pvalue,
                        col_formats['lowCI'][1] % lowCI,
                        col_formats['upCI'][1] % upCI,
                        col_formats['Delay'][1] % delay,
                    ]

                    print(' '.join(data_values), file=resultFile)
                    completed_rows += 1
                    written_since_flush += 1
                    last_triplet_written = (Xi, Yi, Zi)

                    if written_since_flush >= max(1, int(flush_every)):
                        resultFile.flush()
                        written_since_flush = 0

                    if last_triplet_written is not None and completed_rows % max(1, int(checkpoint_every)) == 0:
                        _write_checkpoint(last_triplet_written, completed_rows, completed=False)
                    
                except Exception as e:
                    print("Error during analysis:", file=sys.stderr)
                    traceback.print_exc()
                    continue
    if written_since_flush > 0:
        resultFile.flush()

    if completed_rows == 0:
        print("No valid triplets found for analysis", file=sys.stderr)

    _write_checkpoint(last_triplet_written, completed_rows, completed=True)

def LLApermuPvalue(X, Y, Z, delayLimit, precisionP, LLA_score):
    """Compute permutation-based p-value for LLA score.
    
    Args:
        X, Y, Z: Masked arrays containing time series data
        delayLimit: Maximum time delay to consider
        precisionP: Number of permutations
        LLA_score: Observed LLA score
        
    Returns:
        float: Two-tailed p-value with +1 correction
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

    # Filter out NaN values from permutation results
    PP_set = PP_set[~np.isnan(PP_set)]
    
    # Calculate two-tailed p-value with +1 correction
    # Count how many permuted |scores| >= |observed score|
    b = np.count_nonzero(np.abs(PP_set) >= abs(LLA_score))
    p = (b + 1) / (len(PP_set) + 1)
    
    # Ensure p-value doesn't exceed 1.0
    return float(np.minimum(1.0, p))

def LLAbootstrapCI(X, Y, Z, LLA_score, delayLimit, bootCI, bootNum):
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
        return LLA_score, LLA_score, LLA_score
        
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


def LLAtheoPvalue(LLA_score, series_length, delayLimit, precisionP, x_sd=1.0):
    """Approximate an LLA p-value with the LSA theoretical lookup table.

    The current LLA pipeline is used with delayLimit=0 in our benchmark setup,
    so the score is treated as the range of a standardized partial-sum process
    over the product sequence after normalization.
    """
    precisionP = max(1, int(abs(precisionP)))
    P_table = lsalib.theoPvalue(
        Rmax=series_length,
        Dmax=delayLimit,
        precision=1.0 / float(precisionP),
        x_decimal=lsalib.my_decimal,
    )
    observed_range = abs(LLA_score) * float(series_length)
    return float(
        lsalib.readPvalue(
            P_table,
            R=observed_range,
            N=series_length,
            x_sd=x_sd,
            M=1.0,
            alpha=1.0,
            beta=1.0,
            x_decimal=lsalib.my_decimal,
        )
    )