.. |Logo| image:: https://bitbucket.org/charade/elsa/raw/master/doc/images/elsa_logo.png
   :alt: logo.png
   :height: 50px
   :width: 100px

.. |Pipeline| image:: https://bitbucket.org/charade/elsa/raw/master/doc/images/elsa_pipeline.png
   :alt: lsa_pipeline.png
   :height: 450px
   :width: 540px

.. |elaPipeline| image:: https://bitbucket.org/charade/elsa/raw/master/doc/images/ela_pipeline.png
   :alt: ela_pipeline.png
   :height: 450px
   :width: 540px

|Logo| 

ELSA - Finding Time-Dependent Associations (Now Augmented with Liquid Associaiton)
==========================================================================================

QUICK LINKS
-----------

`Extended Local Similarity Analysis (ELSA) Manuals <https://bitbucket.org/charade/elsa/wiki/Manual>`__

`Extended Liquid Associaiton Analysis (ELA) Manuals <https://bitbucket.org/charade/elsa/wiki/Manual_ela>`__

`FAQ <https://bitbucket.org/charade/elsa/wiki/FAQ>`__

INSTRODUCTION
--------------

In recent years, advances in molecular technologies have empowered researchers with the ability to spatially and temporally characterize natural microbial communities without lab cultivation. Mining and analyzing co-occurrence patterns in these new datasets are fundamental to revealing existing symbiosis relations and microbe-environment interactions. Time series data, in particular, are receiving more and more attention, since not only undirected but also directed associations can be inferred from these datasets.

Researchers typically use techniques like principal component analysis (PCA), multidimensional scaling (MDS), discriminant function analysis (DFA) and canonical correlation analysis (CCA)) to analyze microbial community data under various conditions. Different from these methods, the Local Similarity Analysis (LSA) technique is unique to capture the time-dependent associations (possibly time-shifted) between microbes and between microbe and environmental factors. Significant LSA associations can be interpreted as a partially directed association network for further network-based analysis.

ALso, another limitation of previous investigations was their restriction toward pairwise relationships, while there was no lack of abundance of third-party mediated functionalities in microbial community. Liquid Association analysis was originally developed by K.C. Li for characterizing the internal evolution of expression pattern of a pair of genes (𝑋, 𝑌) depending on a ‘scouting’ (mediator) gene Z. In the newly implemented ELSA, we extended Liquid Association and combined it with LSA to explore the mediated co-varying microbial dynamics, which provides a novel perspective and a useful toolkit to analyze and interpret microbial community time series data. 


METHODS
-------------

|Pipeline|

Figure 1. The analysis workflow of Local Similarity Analysis (LSA) tools. Users start with raw data (matrices of time series) as input and specify their requirements as parameters. The LSA tools subsequently F-transform and normalize the raw data and then calculate the Local Similarity (LS) Scores and the Pearson’s Correlation Coefficients. The tools then assess the statistical significance (P-values) of these correlation statistics using permutation test and filter out insignificant results. Finally, the tools construct a partially directed association network from significant associations.

|elaPipeline|

Figure 2. Liquid Association / Mediated correlation and example Cytoscape diagrams for all liquid association types of factors X, Y and Z: (A) High Z level enhances the positive correlation between X and Y; (B) Low Z level enhances the negative correlation between X and Y; (C) Low Z level enhances the positive correlation between X and Y; (D) High Z level enhances the negative correlation between X and Y. And (E) A flowchart for incorporating Liquid Association (LA) with Local Similarity (LS) Analysis (LSA). First we use LSA to find candidate local and without time-delayed associations between factors X and Y. The results were filtered based on p-values, q-values and effect (LS score). Then, given the significant LSA factors X and Y, we compute LA score to scout any environmental/OTU factors to discover potential mediating factor Z. Next, a permutation test for liquid association is performed and the results were filtered based on p-values, q-values and effect (LA score) to remove insignificant triplets. Finally, we use the software Cytoscape to visualize the results.

SOFTWARE
-------------
    Extended Local Similarity Analysis(ELSA)

    Currently the package is maintained only for Linux (Ubuntu). 
    For other platforms, see the DOCKER section for use.

    More current information of this package is available at:
    http://github.com/labxscut/elsa

DOCKER 
---------------------------------------------

  A Dockerfile is provided to build elsa enabled docker image from a standard Ubuntu docker image. 
  To build a docker image using 'docker build $ELSAPKG', where $ELSAPKG is the unzipped path of elsa.

    ::
      
      host> docker build --no-cache -t labx-elsa .

  Name the built container as your:container; Then mount current data directory to /var/data accessible by docker:

    ::

      host> docker run --name elsa -it -v `pwd`:/var/data/ labx-elsa 
      elsa> lsa_compute --help

INSTALL
-----------------


    (1). # for use #

    First, install the prerequisites: C++ (build-essential), Python(dev), 
    conda, numpy, scipy and matplotlib as specified in setup.py,
    create and activate a virtual env named elsa. 

    Same for Linux, Mac OS, Windows Subsystem Linux, with a Conda-like Virtual Env

    Then, download and unzip the latest master branch to the elsa folder.

    ::

        elsa>  pip install .                            # or python setup.py install
        elsa>  cd test && . test.sh                     # a test script is available


    (2). # for development #

    eLSA is open source and your contributions are greatly welcome.

    First, use git to fork a copy of elsa on github.com:

    ::
        
        elsa> git clone ssh://git@github.com/your/elsa elsa

    Then, make your edits and create a pull request to merge back.

EXECUTABLES
--------------------

    1. The following executables will be available from your python scripts directory (typically already in $PATH).
    2. Use '$script_name -h' for its usage.
    3. A few examples are available in 'test/test.sh' and explained there.

    ::

      lsa_compute                       # main script for LSA/LTA computation to generat the .lsa result
      lsa_version                       # display the git commit hash, i.e. version, of installed esla

      lsa_chkdat                        # provided as is, check the input file's format compatibility
      lsa_query                         # provided as is, query .lsa result and generate network (requires rpy)
      lsa_infer                         # provided as is, plot pairwise .lsa result (requires rpy)
      lsa_sim                           # provided as is, simulate a pair of time series
      lsa_para                          # provided as is, generate parallel jobs for a large input file 
      lsa_fixqv                         # provided as is, fix the q-values for the merged parallel .lsa results 
      lsa_totrend                       # provided as is, convert original series to trend series
    

NOTES
----------------------
    
    The historical R version of LSA is no longer offered, all its computation capacity is avaiblable through eLSA.


CONTACT
----------------------

    lcxia at scut dot edu dot cn

CITATIONS
----------------------

Please cite the references 1, 2, 3, 4 and 5 if the LSA or LTA algorithms and statistical theories was used or reviewed in your work. 
Please cite the references 2 and 3 if any scripts of the ELSA software package was used in your work. 
Please also cite the reference 4 if local trend analysis (LTA) scripts was used in your work. 
Please also cite the reference 6 if local liquid association analysis (LLA) scripts was used in your work. 

    1. (1) Dongmei Ai, Lulu Chen, Jiemin Xie, Longwei Cheng, Fang Zhang, Yihui Luan, Yang Li, Shengwei Hou, Fengzhu Sun, Li Charlie Xia. Identifying local associations in biological time series: algorithms, statistical significance, and applications. Briefings in Bioinformatics 2023, 24(6):bbad390. (https://doi.org/10.1093/bib/bbad390) 
    2. (2) Li C Xia, Dongmei Ai, Jacob Cram, Jed A Fuhrman, Fengzhu Sun. Efficient Statistical Significance Approximation for Local Association Analysis of High-Throughput Time Series Data. Bioinformatics 2013, 29(2):230-237. (https://doi.org/10.1093/bioinformatics/bts668)
    3. (3) Li C Xia, Joshua A Steele, Jacob A Cram, Zoe G Cardon, Sheri L Simmons, Joseph J Vallino, Jed A Fuhrman and Fengzhu Sun. Extended local similarity analysis (eLSA) of microbial community and other time series data with replicates. BMC Systems Biology 2011, 5(S2):S15 (https://doi.org/10.1186/1752-0509-5-S2-S15)
    4. (4) Li C Xia, Dongmei Ai, Jacob Cram, Xiaoyi Liang, Jed Fuhrman, Fengzhu Sun. Statistical significance approximation in local trend analysis of high-throughput time-series data using the theory of Markov chains. BMC Bioinformatics 2015, 16, 301 (https://doi.org/10.1186/s12859-015-0732-8)
    5. (5) Quansong Ruan, Debojyoti Dutta, Michael S. Schwalbach, Joshua A. Steele, Jed A. Fuhrman and Fengzhu Sun Local similarity analysis reveals unique associations among marine bacterioplankton species and environmental factors Bioinformatics 2006, 22(20):2532-2538
    6. (6) Dongmei Ai, Xiaoxin Li, Hongfei Pan, Li Charlie Xia#. Explore Mediated Co-varying Dynamics in Microbial Community using Integrated Local Similarity and Liquid Association Analysis. Accepted bu APBC, to appear in BMC Genomics (2019).
