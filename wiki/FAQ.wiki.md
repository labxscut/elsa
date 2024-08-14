### FAQ

1.  *Range of LS score?*\
    -   As explained in our papers, the LS score itself is not bounded
        between -1 to 1. With in one dataset, one way to compare the
        relative strength of the association is to standardize it by
        dividing all LS scores by the absolute value of the max/min of
        all scores, such that the score become -1 and 1 bounded. The
        statistical significance of the association is however not
        affected. To compare aross datasets is more tricky. One
        potential way is to transform the LS scores into qunatiles with
        signs unchanged, such that it is between -1 and 1. Then the
        strength can be compared in terms of the association\'s relative
        strength amon datasets.
2.  *How the fillMethod works?*\
    -   For fillMethod in eLSA:\
        none: keep na until the data is normalized, then fill it up with
        0\
        others: pass the fillMethod to the interpld function in scipy,\
        see
        \"<http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html>\"
        for how it works.\
3.  *Why tab delimited file produced by Excel not accepted by eLSA
    program?*\
    -   Most likely it is an issue due to different carriage return
        character across platforms. Two methods to go:\
        install dos2unix and convert the input file to unix format\
        issue an \'%s/\^M/\\r/g\' in VI and save as a new file\
    -   Update: \"check\_data\" is now available for checking your input
        format consistency with eLSA expectation.
4.  *The alogrithm has been running about 5 days now - it does appear to
    be runnning smoothly - i submitted the input file i sent you i.e 715
    factors, 39 samples would you suggest cutting this file down?
    shouldl the program run this long?\
    also - how does the program work ? would it go faster if i allocate
    more memory and can it be run on parallel processors? *\
    -   The program is currently running single process. It would be
        helpful if you can run it natively from Linux OS with a faster
        processor. Memory is typically not a concern. How many
        permutations (P) you put? What is the bootstrap number (B). Time
        \~ (B+P)\*F\^2\*N. For 715 factors (F) and 39 samples (N) at
        P=1000 it probably will take a while, depends on your processor
        speed. One way is you can reduce either B, P. Another way is to
        use the bleeding version of the LSA program, which provides
        basic parallel support. You can download the bleeding release
        from the source code site and install it on your blast. Subset
        your input data into, say 10, trunk files and run eLSA on each
        pair (i,j: 0\<i\<=j\<10) of trunk files in the upper triangle
        (so i\<=j) and specify the parameter to tell the program if the
        pair is on diagonal line (i==j). The you can concatenate all the
        result files to get the whole result. This probably can be
        automated by scripting, which I haven\'t got time to write yet.
        One thing is that, you have to recalculated all the Q-values
        after you get all the results. You can collect all the P-values
        and using the R package Qvalue to do it.\
    -   Update: My suggestions are reduce the precision and OTU number
        and try to get a hold of how long it takes to run your analysis
        using the formula provided above before serious jobs. There are
        also \"par\_ana\" and \"fix\_qv\" commands available in the
        bleeding version which divide the jobs and run them parallel if
        you have a PBS system and access to multiple nodes.
5.  *problem with q-values*\
    -   Sometimes the input set of p-values cannot be used to estimate
        q-values successfully using the default storeyQvalue function
        implementation using the interpolation core from scipy. The same
        might also be found using R\'s qvalue package. This is because,
        in some cases, all the pairwise p-values are quite small, a
        setting where the assumption of storeyQvalue and R\'s qvalue -
        the majority p-values is non-significant for q-values so it can
        be used to estimate the uniform rate does not hold. However,
        q-values can be calculated in different ways, while the
        storeyQvalue and R\'s qvalue package all basically implemented
        J.D.Storey\'s method in *Statistical significance for genomewide
        studies PNAS 2003*. There are other methods available but not
        implemented in the eLSA package. A few of them are listed at
        eLSA\'s homepage <http://meta.usc.edu/softs/lsa> . The current
        trade-off is to continue on LSA computation without the q-values
        (putting na in). If needed, p-values can be collected and sent
        to other tools for q-values calculation. The new version is
        pushed to bitbucket and tagged as bleeding.\
    -   Note important if python q-value is not working most likely R\'s
        qvalue package will not work either. The only difference between
        them is the interpolation function used. One is from scipy and
        the other is from R lib.
6.  What is \"-p mix\"? or what is the mixed approach for p-value
    estimation?\
    -   \"mix\" is a combination of theoretical and permutation methods.
        Theoretical p-values are not perfectly linearly correlated to p
        values of permutation methods, thus it is true that \"mix\"
        p-values are not linearly correlated to permutation methods
        either. \"mix\" method saves your time by employing some
        \"intelligence\" - only passing into permutation evaluation of
        Pperm when the Ptheo is significant (\<0.05). The p-values from
        this method are also trusted, and the conclusion you get based
        on the this method should be correct. If you still bother the
        difference between \"mix\" and permutation-only, you can always
        use the permutation-only method.
7.  How to cut off using \"P\" and \"Q\" values?
    -   Cutting off on q-value is basically cutting off on p-value with
        a different strategy and an easier explanation. And if you apply
        both P and Q cutoffs, you get the combined filtering effect:
        **less than probability P to occur by random and with
        probability Q to be false positive**. Basically q values are the
        transformation of p values. So they are telling the same story,
        but with different perspective and interpretations. You can use
        either P values or Q values or a combination of them, whichever
        you think fits your needs best. If you don\'t know what to
        choose, read relevant papers, and find out the conventional
        cut-offs they use and have a try yourself. For more information
        about P and Q values, see Storey et al. *Statistical
        significance for genomewide studies* PNAS 2003.
8.  Have install eror with rpy2 related problem?
    -   rpy2 is no longer a dependency.
    -   If you have previously installed eLSA, completely remove it !!!
        (find it at your python site-packages)
    -   Git clone the latest \"master\" branch or download from
        <https://bitbucket.org/charade/elsa/downloads> (Go for Tags
        version \>v1.0.4)
    -   Install ELSA: *python setup.py install*
    -   Try lsa\_compute see if it works
9.  How many time points are needed to use fast theoretical p-values
    -   you can use theoretical p-values with time points \>10 with no
        delay and \>20 with delay. Please refer to Xia et al. *Efficient
        statistical significance approximation for local similarity
        analysis of high-throughput time series data* Bioinformatics
        29(2):230-237 for details
10. How to choose among normalization schemes?
    -   we recommend the percentile normalization + Z-score
        normalization (i.e. percentileZ option in the command line)
        unless you have other concerns. Percentile is simply rank
        normalization. PercentileZ applies a Z-score transformation
        after Percentile normalization (this is described in recent
        papers). robustZ uses the median instead of mean and mad instead
        of standard deviation for Z-score, otherwise the same as
        PercentileZ. pnz is percentile normalization without zeros.
