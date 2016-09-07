README.rst
------------

INSTRODUCTION
==============

    Extended Local Similarity Analysis(eLSA)
    Currently the package works for Linux (Ubuntu). 
    It might work for Windows with Cygwin (not tested).
    More current information of this package is available @
    http://bitbucket.org/charade/elsa
    
    eLSA Wiki (must read and welcome to contribute) is available @
    http://bitbucket.org/charade/elsa/wiki/Home

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

        sve_vpy> python setup.py install

      (1.5) Now the ELSA executables will be available from "$PWD/sve_vpy/bin". Because you installed ELSA via virtualenv, remember to activate the virtualenv first every time you use ELSA. Also export the environmental variable $ELSA_BIN=$PWD/sve_vpy/bin

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

CONTACT
=============

    lixia at stanford dot edu
