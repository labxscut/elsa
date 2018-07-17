.. |Logo| image:: https://bitbucket.org/charade/elsa/raw/master/doc/images/elsa_logo.png
   :alt: logo.png
   :height: 50px
   :width: 100px

.. |Pipeline| image:: https://bitbucket.org/charade/elsa/raw/master/doc/images/elsa_pipeline.png
   :alt: pipeline.png
   :height: 450px
   :width: 540px

|Logo| ELSA - Finding Time-Dependent Associations in Time Series Datasets 
==========================================================================================

QUICK LINKS
-----------

`Manuals <https://bitbucket.org/charade/elsa/wiki/Manual>`__

`FAQ <https://bitbucket.org/charade/elsa/wiki/FAQ>`__

INSTRODUCTION
==============

    In recent years, advances in molecular technologies have empowered researchers with the ability to spatially and temporally characterize natural microbial communities without lab cultivation (Fuhrman, 2009). Mining and analyzing co-occurrence patterns in these new datasets are fundamental to revealing existing symbiosis relations and microbe-environment interactions (Chaffron et al., 2010; Steele et al., 2011). Time series data, in particular, are receiving more and more attention, since not only undirected but also directed associations can be inferred from these datasets.

    Researchers typically use techniques like principal component analysis (PCA), multidimensional scaling (MDS), discriminant function analysis (DFA) and canonical correlation analysis (CCA)) to analyze microbial community data under various conditions. Different from these methods, the Local Similarity Analysis (LSA) technique is unique to capture the time-dependent associations (possibly time-shifted) between microbes and between microbe and environmental factors (Ruan et al., 2006). Significant LSA associations can be interpreted as a partially directed association network for further network-based analysis.

    Studies adopting the LSA technique have shown interesting and novel discoveries for microbial communities (Paver et al., 2010; Shade et al., 2010; Beman et al., 2011; Steele et al., 2011). However current dataset scale has outdated the old script. To improve computation efficiency, incorporate new features, such as time series data with replicates, and make the analysis technique more accessible to users, we have re-implemented the LSA algorithm as a C++ extension to Python. We also integrated the new LSA tool set with the popular Galaxy framework (Goecks et al., 2010) for web based pipeline analysis.

    Extended Local Similarity Analysis(eLSA)
    Currently the package works for Linux (Ubuntu). 
    It might work for Windows with Cygwin (not tested).
    More current information of this package is available @
    http://bitbucket.org/charade/elsa
    
    eLSA Wiki (must read and welcome to contribute) is available @
    http://bitbucket.org/charade/elsa/wiki/Home

METHODS
==============

|Pipeline|

Figure 1. The analysis workflow of Local Similarity Analysis (LSA) tools. Users start with raw data (matrices of time series) as input and specify their requirements as parameters. The LSA tools subsequently F-transform and normalize the raw data and then calculate the Local Similarity (LS) Scores and the Pearson’s Correlation Coefficients. The tools then assess the statistical significance (P-values) of these correlation statistics using permutation test and filter out insignificant results. Finally, the tools construct a partially directed association network from significant associations. 

DEPENDENCIES
=============

    C++ build environment
        e.g. build-essential and libstdc++6 in Ubuntu
    Python(>=2.7) 
        download @ http://www.python.org/
    Numpy(>=1.0)
        download @ http://www.scipy.org/
    Scipy(>=0.6)
        download @ http://www.scipy.org/
    
    For setting up the dependencies, users may refer to the author's development document @
    http://dl.dropbox.com/u/35182955/Ubuntu_development_environment.html

DOCKER (Platform Independent and Preferred)
=============================================

  A Dockerfile is provided to build elsa enabled docker image from a standard Ubuntu docker image. 
  To build a docker image using 'docker build $ELSAPKG', where $ELSAPKG is the unzipped path of elsa.
  Or download the Dockerfile directly at:

    https://bitbucket.org/charade/elsa/raw/master/Dockerfile

  Name the built container as your:container; Then mount current data directory to /var/data accessible by docker:

    ::

      sudo docker run -it -v `pwd`:/var/data/ your:container
      sudo docker run cd /var/data/ && lsa_compute ...

INSTALL
============


    [Prerequisites]

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

    [Development]

    eLSA is open source and the version controlled repository is @:
        https://bitbucket.org/charade/elsa.
    Use git (http://github.org) to clone a local copy:
        $git clone ssh://git@bitbucket.org/charade/elsa elsa

    Follow standard python module setup to install:
        $cd elsa
        $python setup.py install

    [VirtualBox (Deprecated)]
    The procedure is similar to QIIME VirtualBox install,
        see http://qiime.org/install/virtual_box.html.

    1. Download and install the VirtualBox (VB) version for your machine,
        at http://www.virtualbox.org

    2. Download the SunLab Virtual Box,
        at http://meta.usc.edu/softs/vbox/SunLab.vdi.tgz
        This file is large so it may take
        between a few minutes and a few hours depending on your Internet
  connection speed. You will need to unzip this file, which you can typically do by
        double-clicking on it.

    3. Create a new virtual machine:
        Launch VirtualBox, and create a new machine (press the New button).
        A new window will show up. Click ‘Next’.

        In this screen type SunLab as the name for the virtual machine. Then
        select Linux as the Operating System, and Ubuntu as the version.
        Click Next.

        Select the amount of RAM (memory). You will need at least 512MB, but
        the best option is based on your machine. After selecting the amount of RAM,
        click Next.

        Select “Use existing hard drive”, and click the folder icon next to
        the selector (it has a green up arrow). In the new window click ‘Add’, and
        locate the virtual hard drive that was downloaded in step 2. Click Select and
        then click Next.

        In the new window click Finish.

    4. Double click on the new virtual machine created – it will be called SunLab
        – to boot it for the first time. The default username and password is:
  user

    5. Review any messages that are shown, and select whatever options are best
        for you.

EXECUTABLES
=============

    lsa_compute

USAGE HELP
=============

    (i) Above executables will be available from your python scripts directory.
      Use '-h' to read individual script usage.
    (ii) A simple test example is available at 'test/test.sh' and explained within.

NOTES
=============
    
    A historical R version is available through Prof. Fengzhu Sun's page and is not supported any longer.
    In case the integrated q-value does not work for you, there are many other independent false discovery rate calculation packages, such as locfdr, mixfdr, fuzzyFDR, pi0, fdrci, nFDR.


CONTACT
=============

    fsun at usc dot edu and/or lixia at stanford dot edu

CITATIONS
=============

Please cite the references 1 and 2 if the eLSA python package was used in your study. Please also cite 3 if local trend analysis was used in your study. Please also cite the reference 4 if you used the old R script, which is no loger maintained.

    1. Li C Xia, Dongmei Ai, Jacob Cram, Jed A Fuhrman, Fengzhu Sun. Efficient Statistical Significance Approximation for Local Association Analysis of High-Throughput Time Series Data. Bioinformatics 2013, 29(2):230-237. (https://doi.org/10.1093/bioinformatics/bts668)
    2. Li C Xia, Joshua A Steele, Jacob A Cram, Zoe G Cardon, Sheri L Simmons, Joseph J Vallino, Jed A Fuhrman and Fengzhu Sun. Extended local similarity analysis (eLSA) of microbial community and other time series data with replicates. BMC Systems Biology 2011, 5(S2):S15 (https://doi.org/10.1186/1752-0509-5-S2-S15)
    3. Li C Xia, Dongmei Ai, Jacob Cram, Xiaoyi Liang, Jed Fuhrman, Fengzhu Sun. Statistical significance approximation in local trend analysis of high-throughput time-series data using the theory of Markov chains. BMC Bioinformatics 2015, 16, 301 (https://doi.org/10.1186/s12859-015-0732-8)
    4. Joshua A Steele, Peter D Countway, Li Xia, Patrick D Vigil, J Michael Beman, Diane Y Kim, Cheryl-Emiliane T Chow, Rohan Sachdeva, Adriane C Jones, Michael S Schwalbach, Julie M Rose, Ian Hewson, Anand Patel, Fengzhu Sun, David A Caron, Jed A Fuhrman. Marine bacterial, archaeal and protistan association networks reveal ecological linkages The ISME Journal 2011, 51414–1425
    5. Quansong Ruan, Debojyoti Dutta, Michael S. Schwalbach, Joshua A. Steele, Jed A. Fuhrman and Fengzhu Sun Local similarity analysis reveals unique associations among marine bacterioplankton species and environmental factors Bioinformatics 2006, 22(20):2532-2538
