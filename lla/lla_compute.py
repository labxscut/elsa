#!/usr/bin/env python
#lsa-compute -- computation script for LSA package to perform lsa table calculation 

#License: BSD

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

#public libs
import sys, csv, re, os, time, argparse, string, tempfile, traceback
#numeric libs
import numpy as np
import scipy as sp
sys.path.insert(0, "/workspace")
sys.path.insert(0, "/workspace/lla")
print(">>> sys.path 前三项:", sys.path[:3], file=sys.stderr)

try:
    # When running as installed package
    from lsa import lsalib
    from lsa import compcore
    from lsa.lsalib import noZeroNormalize, percentileNormalize, noneNormalize
    from lsa.lsalib import simpleAverage, sdAverage, simpleMedian, madMedian
    from lsa.lsalib import fillMissing, ma_average, percentileZNormalize
    from lla import llalib
except ImportError:
    # When running in debug mode
    try:
        from . import lsalib
        from .lsalib import noZeroNormalize, percentileNormalize, noneNormalize
        from .lsalib import simpleAverage, sdAverage, simpleMedian, madMedian
        from .lsalib import fillMissing, ma_average
        from . import llalib
    except ImportError:
        import lsalib
        from lsalib import noZeroNormalize, percentileNormalize, noneNormalize
        from lsalib import simpleAverage, sdAverage, simpleMedian, madMedian
        from lsalib import fillMissing, ma_average, percentileZNormalize
        import llalib

print(
    f">>> compcore backend: {getattr(compcore, 'BACKEND_NAME', 'unknown')} "
    f"(using_gpu={getattr(compcore, 'USING_GPU', False)})",
    file=sys.stderr,
)

def validate_input_dimensions(data, repNum, spotNum):
    """Validate input data dimensions and format."""
    if data is None or len(data.shape) != 2:
        raise ValueError(f"Invalid input data shape: {data.shape if data is not None else None}")
        
    expected_cols = spotNum * repNum
    if data.shape[1] != expected_cols:
        raise ValueError(f"Data has {data.shape[1]} columns but expected {expected_cols} "
                        f"(spotNum={spotNum} × repNum={repNum})")


def main():
    start_time = time.time()
    
    __script__ = "lla_compute"
    version_desc = lsalib.safeCmd('lsa_version')
    version_print = "%s (rev: %s) - copyright Li Charlie Xia, lcxia@scut.edu.cn" \
        % (__script__, version_desc)
    print(version_print, file=sys.stderr)
    print("Local Liquid Association Analysis Tool", file=sys.stderr)

    # define arguments: delayLimit, fillMethod, pvalueMethod
    parser = argparse.ArgumentParser(description="Local Liquid Association Analysis Tool")

    parser.add_argument("dataFile", metavar="dataFile", type=str,
                       help="input data file: m by (r * s) tab delimited text; \n \
                             first row starts with '#' as header; \n \
                             m variables, r replicates, s time spots")
    parser.add_argument("resultFile", metavar="resultFile", type=str,
                       help="output result file")
    parser.add_argument("-d", "--delayLimit", dest="delayLimit", default=3, type=int,
                       help="maximum time delay (default: 3, range: 0-6)")
    parser.add_argument("-p", "--pvalueMethod", dest="pvalueMethod", default="perm",
                       help="p-value calculation method (default: perm; supported: perm, theo, mix)")
    parser.add_argument("-x", "--precision", dest="precision", default=1000, type=int,
                       help="precision for p-value calculation (default: 1000)")
    parser.add_argument("-b", "--bootNum", dest="bootNum", default=0, type=int,
                       choices=[0, 100, 200, 500, 1000, 2000],
                       help="bootstrap iterations for CI (default: 0)")
    parser.add_argument("-r", "--repNum", dest="repNum", default=1, type=int,
                       help="replicates per time spot (default: 1)")
    parser.add_argument("-s", "--spotNum", dest="spotNum", default=4, type=int,
                       help="number of time spots (default: 4)")
    parser.add_argument("-m", "--minOccur", dest="minOccur", default=0.50, type=float,
                       help="minimum occurrence threshold (default: 0.50)")
    parser.add_argument("-t", "--transFunc", dest="transFunc", default='simple',
                       choices=['simple', 'SD', 'Med', 'MAD'],
                       help="replicate summarization method (default: simple)")
    parser.add_argument("-k", "--keep-trace", dest="keep_trace", action='store_true',
                    help="Enable trace output for LLA analysis (records optimal path indices; memory intensive)")
    parser.add_argument("-f", "--fillMethod", dest="fillMethod", default='linear',
                       choices=['none', 'zero', 'linear', 'quadratic', 'cubic', 'slinear', 'nearest'],
                       help="missing value fill method (default: linear)")
    parser.add_argument("-n", "--normMethod", dest="normMethod", default='pnz',
                       choices=['percentile', 'pnz', 'none'],
                       help="data normalization method (default: pnz)")
    parser.add_argument("--resume", dest="resume", action='store_true',
                        help="resume from checkpoint and append to existing result file")
    parser.add_argument("--checkpoint-file", dest="checkpoint_file", default=None,
                        help="checkpoint path (default: <resultFile>.ckpt.json)")
    parser.add_argument("--flush-every", dest="flush_every", default=100, type=int,
                        help="flush output file every N rows (default: 100)")
    parser.add_argument("--checkpoint-every", dest="checkpoint_every", default=1000, type=int,
                        help="save checkpoint every N rows (default: 1000)")

    args = parser.parse_args()

    # Print parameters
    checkpoint_file = args.checkpoint_file or (args.resultFile + ".ckpt.json")

    print("\t".join(['delayLimit','fillMethod','pvalueMethod','dataFile','resultFile',
                    'repNum','spotNum','bootNum','transFunc','normMethod','precision','keep_trace',
                    'resume','checkpoint_file','flush_every','checkpoint_every']))
    print("\t".join(['%s']*16) % (args.delayLimit, args.fillMethod, args.pvalueMethod,
                                 args.dataFile, args.resultFile, args.repNum,
                                 args.spotNum, args.bootNum, args.transFunc,
                                 args.normMethod, args.precision, args.keep_trace,
                                 args.resume, checkpoint_file, args.flush_every, args.checkpoint_every))

    try:
        data_handle = open(args.dataFile, 'r')
        result_mode = 'a' if args.resume else 'w'
        result_handle = open(args.resultFile, result_mode)

        # Map transformation function names to actual functions
        transform_funcs = {
            'simple': lsalib.simpleAverage,
            'SD': lsalib.sdAverage,
            'Med': lsalib.simpleMedian,
            'MAD': lsalib.madMedian
        }
        
        normalize_funcs = {
            'none': lsalib.noneNormalize,
            'percentile': lsalib.percentileNormalize,
            'pnz': lsalib.percentileZNormalize
        }
        
        # Get the appropriate functions based on args
        fTransform = transform_funcs.get(args.transFunc)
        if fTransform is None:
            raise ValueError(f"Invalid transform function: {args.transFunc}")
            
        zNormalize = normalize_funcs.get(args.normMethod)
        if zNormalize is None:
            raise ValueError(f"Invalid normalization method: {args.normMethod}")

        # Read and validate input data
        firstData = np.genfromtxt(data_handle, comments='#', delimiter='\t',
                               missing_values=['na','','NA'],
                               filling_values=np.nan,
                               usecols=list(range(1, args.spotNum*args.repNum+1)))
        
        validate_input_dimensions(firstData, args.repNum, args.spotNum)
        
        # Read factor labels
        data_handle.seek(0)
        factorLabels = np.genfromtxt(data_handle, comments='#', delimiter='\t',
                                  usecols=[0], dtype=str).tolist()
        factorNum = firstData.shape[0]
        
        # Create masked array and reshape
        cleanData = np.ma.zeros((factorNum, args.repNum, args.spotNum), dtype='float')
        for i in range(factorNum):
            for j in range(args.repNum):
                series = firstData[i][j::args.repNum]
                filled = fillMissing(series, args.fillMethod)
                cleanData[i,j] = np.ma.array(filled, mask=np.isnan(filled))
        # cleanData[i, j, k] = ith factor, jth replicates, kth time spots.
        
        # Run analysis with transformed data and specified functions
        llalib.applyLLAnalysis(cleanData, factorLabels,
                            delayLimit=args.delayLimit,
                            bootCI=0.95,
                            bootNum=args.bootNum,
                            minOccur=args.minOccur,
                            pvalueMethod=args.pvalueMethod,
                            precision=args.precision,
                            fillMethod=args.fillMethod,
                            normMethod=args.normMethod,
                            fTransform=fTransform,
                            zNormalize=zNormalize,
                            resultFile=result_handle,
                            keep_trace=args.keep_trace,
                            checkpoint_file=checkpoint_file,
                            resume=args.resume,
                            flush_every=args.flush_every,
                            checkpoint_every=args.checkpoint_every)
                            
    except Exception as e:
        print("Error during analysis:", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)
    finally:
        try:
            data_handle.close()
        except Exception:
            pass
        try:
            result_handle.close()
        except Exception:
            pass

    print("Finishing up...", file=sys.stderr)
    end_time = time.time()
    print(f"Time elapsed {end_time-start_time:.2f} seconds", file=sys.stderr)

if __name__=="__main__":
    main()
