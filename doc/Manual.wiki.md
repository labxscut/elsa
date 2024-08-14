ELSA Manual
-----------

### Input Format (check\_data)

Please use check\_data first to check if your data file is compatible
with ELSA. Transferring text files among Mac, Linux and Windows can
easily mess up your formats. (Here is the reason:
<https://en.wikipedia.org/wiki/Newline>). So always do this before
lsa\_compute.

    usage: check_data [-h] dataFile repNum spotNum

    Auxillary tool to new LSA package for checking data format

    positional arguments:
      dataFile    the data file
      repNum      replicates number
      spotNum     timepoints number

    optional arguments:
      -h, --help  show this help message and exit

The input has to be a tab delimited matrix file, for example the
following one:

    #F3T4R2  t1r1  t1r2  t2r1  t2r2  t3r1  t3r2  t4r1  t4r2
    f1       na    2     3     0     na    1     3     5
    f2       10    na    na    3     na    9     3     3
    f3       -2    -4    na    1     na    0     1     1

For this example file, spotNum=4 and repNum=2.

So each column is one replicate from one time point. t1r1 is replicate
one from timepoint one. And each row is a factor. f1 is factor one.

Make the top left cell whatever but start with an \'\#\'. \'na\' is
reserved for missing value. You might what to take a note of the number
of factors, timespots and replicates, some which are needed for
executing the program.

If you are using Excel for preparing the input file, remember to take
out any trailing and leading empty rows, columns or cells. Make the
table a real rectangle, not only a visually one! That shall do the
input.

### Computation (lsa\_compute)


    lsa_compute (rev: v1.0.2@GIT: fd167ef) - copyright Li Charlie Xia, lixia@stanford.edu
    usage: lsa_compute [-h] [-e EXTRAFILE] [-d DELAYLIMIT] [-m MINOCCUR]
                       [-p {perm,theo,mix}] [-x PRECISION]
                       [-b {0,100,200,500,1000,2000}] [-r REPNUM] [-s SPOTNUM]
                       [-t {simple,SD,Med,MAD}]
                       [-f {none,zero,linear,quadratic,cubic,slinear,nearest}]
                       [-n {percentile,percentileZ,pnz,robustZ,rnz,none}]
                       [-q {scipy}] [-T TRENDTHRESH] [-a APPROXVAR]
                       [-v PROGRESSIVE]
                       dataFile resultFile

    positional arguments:
      dataFile              the input data file, m by (r * s)tab delimited text;
                            top left cell start with '#' to mark this is the
                            header line; m is number of variables, r is number of
                            replicates, s it number of time spots; first row:
                            #header s1r1 s1r2 s2r1 s2r2; second row: x ?.?? ?.??
                            ?.?? ?.??; for a 1 by (2*2) data
      resultFile            the output result file

    optional arguments:
      -h, --help            show this help message and exit
      -e EXTRAFILE, --extraFile EXTRAFILE
                            specify an extra datafile, otherwise the first
                            datafile will be used and only lower triangle entries
                            of pairwise matrix will be computed
      -d DELAYLIMIT, --delayLimit DELAYLIMIT
                            specify the maximum delay possible, default: 0, must
                            be an integer >=0 and <spotNum
      -m MINOCCUR, --minOccur MINOCCUR
                            specify the minimum occurence percentile of all times,
                            default: 50,
      -p {perm,theo,mix}, --pvalueMethod {perm,theo,mix}
                            specify the method for p-value estimation, default:
                            pvalueMethod=perm, i.e. use permutation theo:
                            theoretical approximaton; if used also set -a value.
                            mix: use theoretical approximation for pre-screening
                            if promising (<0.05) then use permutation.
      -x PRECISION, --precision PRECISION
                            permutation/precision, specify the permutation number
                            or precision=1/permutation for p-value estimation.
                            default is 1000, must be an integer >0
      -b {0,100,200,500,1000,2000}, --bootNum {0,100,200,500,1000,2000}
                            specify the number of bootstraps for 95% confidence
                            interval estimation, default: 100, choices: 0, 100,
                            200, 500, 1000, 2000. Setting bootNum=0 avoids
                            bootstrap. Bootstrap is not suitable for non-
                            replicated data.
      -r REPNUM, --repNum REPNUM
                            specify the number of replicates each time spot,
                            default: 1, must be provided and valid.
      -s SPOTNUM, --spotNum SPOTNUM
                            specify the number of time spots, default: 4, must be
                            provided and valid.
      -t {simple,SD,Med,MAD}, --transFunc {simple,SD,Med,MAD}
                            specify the method to summarize replicates data,
                            default: simple, choices: simple, SD, Med, MAD NOTE:
                            simple: simple averaging SD: standard deviation
                            weighted averaging Med: simple Median MAD: median
                            absolute deviation weighted median;
      -f {none,zero,linear,quadratic,cubic,slinear,nearest}, --fillMethod {none,zero,linear,quadratic,cubic,slinear,nearest}
                            specify the method to fill missing, default: none,
                            choices: none, zero, linear, quadratic, cubic,
                            slinear, nearest operation AFTER normalization: none:
                            fill up with zeros ; operation BEFORE normalization:
                            zero: fill up with zero order splines; linear: fill up
                            with linear splines; slinear: fill up with slinear;
                            quadratic: fill up with quadratic spline; cubic: fill
                            up with cubic spline; nearest: fill up with nearest
                            neighbor
      -n {percentile,percentileZ,pnz,robustZ,rnz,none}, --normMethod {percentile,percentileZ,pnz,robustZ,rnz,none}
                            must specify the method to normalize data, default:
                            robustZ, choices: percentile, none, pnz, percentileZ,
                            robustZ or a float NOTE: percentile: percentile
                            normalization, including zeros (only with perm) pnz:
                            percentile normalization, excluding zeros (only with
                            perm) percentileZ: percentile normalization +
                            Z-normalization rnz: percentileZ normalization +
                            excluding zeros + robust estimates (theo, mix, perm
                            OK) robustZ: percentileZ normalization + robust
                            estimates (with perm, mix and theo, and must use this
                            for theo and mix, default)
      -q {scipy}, --qvalueMethod {scipy}
                            specify the qvalue calculation method, scipy: use
                            scipy and storeyQvalue function, default
      -T TRENDTHRESH, --trendThresh TRENDTHRESH
                            if trend series based analysis is desired, use this
                            option NOTE: when this is used, must also supply
                            reasonble values for -p, -a, -n options
      -a APPROXVAR, --approxVar APPROXVAR
                            if use -p theo and -T, must set this value
                            appropriately, precalculated -a {1.25, 0.93, 0.56,0.13
                            } for i.i.d. standard normal null and -T {0, 0.5, 1,
                            2} respectively. For other distribution and -T values,
                            see FAQ and Xia et al. 2013 in reference
      -v PROGRESSIVE, --progressive PROGRESSIVE
                            specify the number of progressive output to save
                            memory, default: 0, 2G memory is required for 1M
                            pairwise comparison.

So we can analyze the above example file by:

      lsa_compute ../test/testna.txt ../test/testna.lsa -r 2 -s 4 -d 1

eLSA will take ../test/testna.txt as input, and knows it has 4 timespots
each with 2 replicates. And eLSA will analyze it with maximum delay of 1
time unit. The output file is explained below.

### Output

    X   Y   LS  lowCI   upCI    Xs  Ys  Len Delay   P   PCC Ppcc    SPCC    Pspcc   Dspcc   SCC Pscc    SSCC    Psscc   Dsscc   Q   Qpcc    Qspcc   Qscc    Qsscc   Xi  Yi
    f1  f2  -0.349677   -0.349677   -0.349677   3   3   2   0   0.316000    -0.512027   0.487973    0.520401    0.651565    1   -0.210819   0.789181    0.500000    0.666667    1   0.451429    0.731960    0.782394    1.000000    1.000000    1   2
    f1  f3  0.349677    0.349677    0.349677    3   3   2   0   0.675000    0.217606    0.782394    0.217606    0.782394    0   0.210819    0.789181    -0.500000   0.666667    1   0.642857    0.782394    0.782394    1.000000    1.000000    1   3
    f2  f3  -1.125007   -1.125007   -1.125007   1   1   4   0   0.215000    -0.827992   0.172008    0.991241    0.084323    -1  -1.000000   nan -1.000000   nan 0   0.451429    0.516024    0.252970    nan nan 2   3

-   X: factor name X
-   Y: factor name Y
-   LS: Local Similarity Score
-   low/upCI: low or up 95% CI for LS
-   Xs: align starts position in X
-   Ys: align starts position in Y
-   Len: align length
-   Delay: calculated delay for align, Xs-Ys
-   P,Q: p/q-value for LS
-   PCC,Ppcc,Qpcc: Pearson\'s Correlation Coefficient, p/q-value for PCC
-   SCC,Pscc,Qscc: Spearman\'s Correlation Coefficient, p/q-value for
    SCC
-   SPCC,Pspcc,Qspcc,Dspcc: delay-Shifted Pearson\'s Correlation
    Coefficient, p/q-value, delay size for SPCC
-   SSCC,Psscc,Qsscc,Dsscc: delay-Shifted Spearman\'s Correlation
    Coefficient, p/q-value, delay size for SSCC

### Speed Up (par\_ana)

You can use par\_ana.py and ssa.py to speed up your analysis using
parallelism in high performance computing clusters.

Then \"par\_ana -h\" tells you how to use the script for computing. In
the singleCmd options, with your normal single line lsa\_comput command,
now replace your input and output by %s symbol. The input and output is
now supplied to multiInput and multiOutput options now. Here the input
is ARISA.txt and the output is ARISA.lsa.

    Example: par_ana ARISA20.txt ARISA20.lsa 'lsa_compute %s %s -e ARISA20.txt -s 127 -r 1 -p theo' $PWD
    Example: par_ana ARISA20.txt ARISA20.la 'la_compute %s ARISA20.laq %s -s 127 -r 1 -p 1000' $PWD
    vmem= 2000mb
    usage: par_ana [-h] [-d DRYRUN] multiInput multiOutput singleCmd workDir

    Multiline Input Split and Combine Tool for LSA and LA

    positional arguments:
      multiInput            the multiline input file
      multiOutput           the multiline output file
      singleCmd             single line command line in quotes
      workDir               set current working directory

    optional arguments:
      -h, --help            show this help message and exit
      -d DRYRUN, --dryRun DRYRUN
                            generate pbs only

par\_ana will use ssa.py to submit the pbs jobs to batch system.

    usage: ssa.py [-h] pbsFile

    MCB Queue Checking and Submission Tool

    positional arguments:
      pbsFile     single pbs file to be submitted

    optional arguments:
      -h, --help  show this help message and exit

Put the ssa.py (shipped in elsa\_pkg/lsa/ssa.py) into your path and set
the queue parameters correctly set inside the script.

Example: you have 63 cores with \#300Gig\# mem in the queue \#main\# and
your username is \#user\#.

    core_max=63
    mem_max=300
    uname="user"
    qname="main"

### FAQ

Wondering: 1. whether to permutation or theoretical p-values? 2. which
normalization to choose?

Or any other doubts, first refer to the [FAQ](FAQ.wiki).

Have fun!
