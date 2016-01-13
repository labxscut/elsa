README.txt

INSTRODUCTION:
    Extended Local Similarity Analysis(eLSA)
    Currently the package works for Linux (Ubuntu) and Mac (Snow Leopard) platforms.
    It might work for Windows with Cygwin (not tested).
    More current information of this package is available @
		http://meta.usc.edu/softs/lsa
    eLSA Wiki (must read and welcome to contribute) is available @
		http://bitbucket.org/charade/elsa/wiki/Home

DEPENDENCIES:
    C++ build environment
        e.g. build-essential and libstdc++6 in Ubuntu and Xcode in Mac 
    Python(>=2.7) 
        download @ http://www.python.org/
    Numpy(>=1.0)
        download @ http://www.scipy.org/
    Scipy(>=0.6)
        download @ http://www.scipy.org/
    Mercurial(optional, required for displaying version number in commandline tools)
	introduction @ http://en.wikipedia.org/wiki/Mercurial
    
    For setting up the dependencies, users may refer to the author's development document @
		http://dl.dropbox.com/u/35182955/Mac_development_environment.html
		http://dl.dropbox.com/u/35182955/Ubuntu_development_environment.html

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
        $tar -zxvf charade-elsa-release.tar.gz
        $cd charade-elsa-release
        $python setup.py install
        $cd test      #test the scripts are workable
        $. test.sh    #ad hoc test of the script on test data

    [MAC] (Snow Leopard)

    Download the latest release of eLSA from http://meta.usc.edu/softs/lsa.
    Follow standard python module setup to install:
        $tar -zxvf charade-elsa-release.tar.gz
        $cd charade-elsa-release
        $python setup.py install

    [Development]

    eLSA is open source and the version controlled repository is @:
    	  https://bitbucket.org/charade/elsa.
    Use mercurial tools (http://mercurial.selenic.com) to download a local copy:
        $hg clone ssh://hg@bitbucket.org/charade/elsa elsa-tip

    Follow standard python module setup to install:
        $cd elsa-tip
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

EXECUTABLES:
    lsa_compute
    lsa_query

USAGE HELP:
    (i) Above executables will be available from your python scripts directory.
    	Use '-h' to read individual script usage.
    (ii) A simple test example is available at 'lsa/test.sh' and explained within.

CONTACT:
    lxia at usc dot edu
