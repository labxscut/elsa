README.txt

INSTRODUCTION:
    Extended Local Similarity Analysis(eLSA)
    Currently the package works for Linux (Ubuntu) and Mac (Snow Leopard) platforms.
    It might work for Windows with Cygwin (not tested).
    More current information of this package is available @
		http://meta.usc.edu/softs/lsa

DEPENDENCIES:
	C++ build environment
		e.g. build-essentials and libstdc++ in Ubuntu and Xcode in Mac 
    Python(>=2.7) 
        download @ http://www.python.org/
    Numpy(>=1.0)
        download @ http://www.scipy.org/
    Scipy(>=0.6)
        download @ http://www.scipy.org/
    
    For setting up the dependencies, users may refer to the author's development document @
		http://dl.dropbox.com/u/35182955/Mac_development_environment.html
		http://dl.dropbox.com/u/35182955/Ubuntu_developement_environment.html

FILES:
    LICENSE.txt
    README.txt
    MANIFEST.in
    setup.py
    lsa/*
    test/*

INSTALL:
	
	[Prerequisites]

    Please fullfill the prerequisites of C++, Python (with development and setuptools),
    numpy, scipy and biopython as described in README.txt before installing eLSA.
    
    [Linux] (Ubuntu)

    Download the latest release of eLSA from http://meta.usc.edu/softs/lsa.
    Follow standard python module setup to install:
        $tar -zxvf charade-lsa-release.tar.gz
        $cd charade-lsa-release
        $python setup.py install

    [MAC] (Snow Leopard)

    Download the latest release of GRAMMy from http://meta.usc.edu/softs/lsa.
    Follow standard python module setup to install:
        $tar -zxvf charade-lsa-release.tar.gz
        $cd charade-lsa-release
        $python setup.py install

    [Development]

    eLSA is open source and the version controlled repository is @:
    	https://bitbucket.org/charade/lsa.
    Use mercurial tools (http://mercurial.selenic.com) to download a local copy:
        $hg clone ssh://hg@bitbucket.org/charade/lsa lsa-tip

    Follow standard python module setup to install:
        $cd lsa-tip
        $python setup.py install

EXECUTABLES:
    lsa_compute
    lsa_query

USAGE HELP:
    (i) Above executables will be available from your python scripts directory.
    	Use '-h' to read individual script usage.
    (ii) A simple test example is available at 'lsa/test.sh' and explained within.

CONTACT:
    lxia@usc.edu
