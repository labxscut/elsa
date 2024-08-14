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
    Currently the package was developed for Linux (Ubuntu). 
    For other platforms, see the DOCKER section for help.

    More current information of this package is available @
    http://bitbucket.org/charade/elsa
    
    ELSA Wiki (a must read and welcome to contribute) is available @
    http://bitbucket.org/charade/elsa/wiki/Home

    Use of external resources:

    C++ build environment
        e.g. build-essential and libstdc++6 in Ubuntu
    Python(>=2.7) 
        download @ http://www.python.org/
    Numpy(>=1.0)
        download @ http://www.scipy.org/
    Scipy(>=0.6)
        download @ http://www.scipy.org/

DOCKER (Platform Independent and Preferred)
---------------------------------------------

  A Dockerfile is provided to build elsa enabled docker image from a standard Ubuntu docker image. 
  To build a docker image using 'docker build $ELSAPKG', where $ELSAPKG is the unzipped path of elsa.
  Or download the Dockerfile directly at https://bitbucket.org/charade/elsa/raw/master/Dockerfile:

  ::
      
      git clone https://bitbucket.org/charade/elsa
      cd elsa && docker build . -t elsa_container

  This built and named the elsa container as elsa_container; 
  Then mount current data directory to /var/data accessible by docker:

  ::
      
      docker run -n elsa1 -it -v `pwd`:/var/data/ elsa_container
      docker exec elsa1 sh -c "cd /var/data/ && lsa_compute ..."


INSTALL
-----------------


    1. Install Prerequisites

    Please fullfill the prerequisites of C++, Python (with development and setuptools),
    numpy, scipy and biopython as described in README.txt before installing eLSA.
    
    [Linux] (e.g. Ubuntu)

    Download the latest master branch of eLSA from https://bitbucket.org/charade/elsa/get/master.tar.gz .
    Follow standard python module setup to install:
        $tar -zxvf charade-elsa-master.tar.gz
        $cd charade-elsa-$your_master_commit_id
        $python setup.py install
        $cd test      #test the scripts are workable
        $. test.sh    #ad hoc test of the script on test data

    [Linux] (virtualenv)

    Install ELSA through system/site python and virtualenv

      **This is the MOST RECOMMENDED WAY for installation**

     (1.1) virtualenv command is standard with Python 2.7 or later. If it is not present, please see https://virtualenv.pypa.io for details to install virtualenv for your python. Possibly as simple as:

      ::

        sudo easy_install pip
        sudo pip install virtualenv

      Ask your IT manager to help install it for you if you have permission difficulties.

      (1.2) When your system python has virtualenv, make sure your $PYTHONPATH is set to empty and follow steps below:

      ::

        >virtualenv-2.7 vpy27 --no-site-packages

      (1.3) Then you can activate this virtual python:

      ::

        >source vpy27/bin/activate
        >pip install numpy
        >pip install scipy

      (1.4) Now under your virtualenv, the dependencies will be automatically setup:

      ::

        vpy27> python setup.py install

      (1.5) Now the ELSA executables will be available from "$PWD/vpy27/bin". Because you installed ELSA via virtualenv, remember to activate the virtualenv first every time you use ELSA. Also export the environmental variable $ELSA_BIN=$PWD/vpy27/bin

    2. Development

    eLSA is open source and the version controlled repository is @:
        https://bitbucket.org/charade/elsa.
    Use git (http://github.org) to clone a local copy:
        $git clone ssh://git@bitbucket.org/charade/elsa elsa

    Follow standard python module setup to install:
        $cd elsa
        $python setup.py install

EXECUTABLES
--------------------

  ::

    lsa_compute                       # for LSA/LTA computation
    la_compute                        # for LA computation

USAGE
---------------------

    (i) Above executables will be available from your python scripts directory.
      Use '-h' to read individual script usage.
    (ii) A simple test example is available at 'test/test.sh' and explained within.

NOTES
----------------------
    
    A historical R version is available through Prof. Fengzhu Sun's page and is not supported any longer.
    In case the integrated q-value does not work for you, there are many other independent false discovery rate calculation packages, such as locfdr, mixfdr, fuzzyFDR, pi0, fdrci, nFDR.


CONTACT
----------------------

    fsun at usc dot edu and/or lixia at stanford dot edu

CITATIONS
----------------------

Please cite the references 1 and 2 if any part of the ELSA python package was used in your study. Please also cite 3 if local trend analysis (LTA) was used in your study. Please also cite reference 6 if extended liquid association analysis (ELA) was used in your study. Please also cite the reference 4 and 5 if you used the old LSA R script, which is no loger maintained. 

    1. Li C Xia, Dongmei Ai, Jacob Cram, Jed A Fuhrman, Fengzhu Sun. Efficient Statistical Significance Approximation for Local Association Analysis of High-Throughput Time Series Data. Bioinformatics 2013, 29(2):230-237. (https://doi.org/10.1093/bioinformatics/bts668)
    2. Li C Xia, Joshua A Steele, Jacob A Cram, Zoe G Cardon, Sheri L Simmons, Joseph J Vallino, Jed A Fuhrman and Fengzhu Sun. Extended local similarity analysis (eLSA) of microbial community and other time series data with replicates. BMC Systems Biology 2011, 5(S2):S15 (https://doi.org/10.1186/1752-0509-5-S2-S15)
    3. Li C Xia, Dongmei Ai, Jacob Cram, Xiaoyi Liang, Jed Fuhrman, Fengzhu Sun. Statistical significance approximation in local trend analysis of high-throughput time-series data using the theory of Markov chains. BMC Bioinformatics 2015, 16, 301 (https://doi.org/10.1186/s12859-015-0732-8)
    4. Joshua A Steele, Peter D Countway, Li Xia, Patrick D Vigil, J Michael Beman, Diane Y Kim, Cheryl-Emiliane T Chow, Rohan Sachdeva, Adriane C Jones, Michael S Schwalbach, Julie M Rose, Ian Hewson, Anand Patel, Fengzhu Sun, David A Caron, Jed A Fuhrman. Marine bacterial, archaeal and protistan association networks reveal ecological linkages The ISME Journal 2011, 51414–1425
    5. Quansong Ruan, Debojyoti Dutta, Michael S. Schwalbach, Joshua A. Steele, Jed A. Fuhrman and Fengzhu Sun Local similarity analysis reveals unique associations among marine bacterioplankton species and environmental factors Bioinformatics 2006, 22(20):2532-2538
    6. Dongmei Ai, Xiaoxin Li, Hongfei Pan, Li Charlie Xia*. Extending Liquid Association to Explore Mediated Co- varying Dynamics in Marine Microbial Community. Manuscript under review (2018).