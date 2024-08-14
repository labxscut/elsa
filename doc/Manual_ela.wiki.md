ELA Manual
----------

### Input Format (check\_data)

Please use check\_data first to check if your data file is compatible
with ELA. Transferring plain text files among Mac, Linux and Windows can
easily mess up your formats. (Here is the reason:
<https://en.wikipedia.org/wiki/Newline>). So always do this before
la\_compute.

    usage: check_data [-h] dataFile repNum spotNum

    Auxillary tool to new ELA package for checking data format

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

### Computation (la\_compute)

    la_compute (rev: v1.0.2) - copyright Li Charlie Xia, lixia@stanford.edu
    usage: la_compute [-h] [-xi XICOL] [-yi YICOL] [-x PRECISION] [-p {perm}]
                      [-m MINOCCUR] [-b {0,100,200,500,1000,2000}] [-r REPNUM]
                      [-s SPOTNUM] [-t {simple,SD,Med,MAD}]
                      [-f {none,zero,linear,quadratic,cubic,slinear,nearest}]
                      [-n NORMMETHOD]
                      dataFile scoutFile resultFile

    Extended Liquid Association Analysis Tools

    positional arguments:
      dataFile              the input data file, m by (r * s)tab delimited text;
                            top left cell start with '#' to mark this is the
                            header line; m is number of variables, r is number of
                            replicates, s it number of time spots; first row:
                            #header s1r1 s1r2 s2r1 s2r2; second row: x ?.?? ?.??
                            ?.?? ?.??; for a 1 by (2*2) data
      scoutFile             the input datafile specify the scouting pairs, it can
                            be any tab delimited file (e.g. .lsa) with (xi, yi)
                            pair indecies for scouting pairs
      resultFile            the output result file

    optional arguments:
      -h, --help            show this help message and exit
      -xi XICOL, --xiCol XICOL
                            specify the x-th column to store Xi indecies
      -yi YICOL, --yiCol YICOL
                            specify the y-th column to store Yi indecies
      -x PRECISION, --precision PRECISION
                            permutation/precision, specify the permutation number
                            or precision=1/permutation for p-value estimation.
                            must be integer >0
      -p {perm}, --pvalueMethod {perm}
                            specify the method for p-value estimation, default:
                            pvalueMethod=perm, i.e. use permutation. 
                            it is the only option available for ELA.
      -m MINOCCUR, --minOccur MINOCCUR
                            specify the minimum occurence percentile of all times,
                            default: 50,
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
      -n NORMMETHOD, --normMethod NORMMETHOD
                            specify the method to normalize data, default:
                            percentile, choices: percentile, none, pnz,
                            percentileZ NOTE: percentile: percentile
                            normalization, including zeros pnz: percentile
                            normalization, excluding zeros percentileZ: percentile
                            normalization + Z-normalization none or a float number
                            for variance: no normalization and calculate Ptheo
                            with user specified variance, default=1

      lsa_compute ../test/testna.txt ../test/testna.lsa -r 2 -s 4 -d 0
      la_compute ../test/testna.txt ../test/testna.lsa ../test/testna.la -r 2 -s 4

In the first step, ELSA will take ../test/testna.txt as input, and knows
it has 4 timespots each with 2 replicates. And eLSA will analyze it with
maximum delay of 0 time unit.

In the next step, ELA will take ../test/testna.txt and
../test/testna.lsa as input, and knows it has 4 timespots each with 2
replicates.

The output file ../test/testna.la is explained below.

### Output

    X     Y   Z   LA          lowCI     upCI        P           Q       Xi  Yi  Zi
    f1  f2  f3  -0.795101   -0.795101   -0.795101   0.320000    0.00    1     2   3

-   X: factor name of X
-   Y: factor name of Y
-   Z: factor name of Z
-   LA: Liquid Association Score
-   low/upCI: low or up 95% CI for LA score
-   P,Q: p/q-value for LA
-   Xi: of X
-   Yi: factor name of Y
-   Zi: factor name of Z

### Speed Up (par\_ana)

Note: these instructions are provisional and not supported.

You can use par\_ana.py and ssa.py to speed up your analysis using
parallelism in high performance computing clusters.

Then \"par\_ana -h\" tells you how to use the script for computing. In
the singleCmd options, with your normal single line lsa\_comput command,
now replace your input and output by %s symbol. The input and output is
now supplied to multiInput and multiOutput options now. Here the input
is ARISA.txt and the output is ARISA.lsa.

    Example: par_ana ARISA20.txt ARISA20.lsa 'la_compute %s %s -e ARISA20.txt -s 127 -r 1 -p theo' $PWD
    Example: par_ana ARISA20.txt ARISA20.la 'la_compute %s ARISA20.laq %s -s 127 -r 1 -p 1000' $PWD
    vmem= 2000mb
    usage: par_ana [-h] [-d DRYRUN] multiInput multiOutput singleCmd workDir

    Multiline Input Split and Combine Tool for ELSA and ELA

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

Have fun!
