README.txt(v0.1.0)


PACKAGE NAME:
    Local Similarity Analysis(LSA)

DEPENDENCY:
    Python(>=2.7) 
        download @ http://www.python.org/
    Numpy(>=1.0)
        download @ http://www.scipy.org/
    Scipy(>=0.6)
        download @ http://www.scipy.org/

FILES:
    LICENSE.txt
    README.txt
    INSTALL.txt
    MANIFEST.in
    setup.py
    lsa/__init__.py
    lsa/lsalib.py
    lsa/lsaio.py
    lsa/compcore.h
    lsa/compcore.cpp
    lsa/compcore.py
    lsa/lsa-compute.py
    lsa/lsa-query.py
    doc/*
    test/*

INSTALL:
    For Linux/Mac/Unix, first install all prequisites, then after untar: 
        python setup.py install
    Note: 
        for details look into INSTALL.txt

EXECUTABLES:
    Unix/Linux: 
        lsa-compute.py and lsa-query.py can be found in "$HOME/bin"
        make sure "$HOME/bin" is in your PATH environment variable

USAGE HELP:
    lsa_compute -h
    lsa_query -h

DOCUMENTATION:
    (i) developement information of lsa package and lsa.lsaio and lsa.lsalib models 
    can be found by issuing "pydoc lsa", "pydoc lsa.lsaio", and "pydoc lsa.lsalib" 
    commands on command line, respectively.
    (ii) information of how to use this package for analysis data by examples can
    be found in documentation in "doc" subdirectory

CONTACT:
    lxia@usc.edu
