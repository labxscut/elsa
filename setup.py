#!/usr/bin/env python

#License: BSD

#Copyright (c) 2008 Li Charles Xia
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions
#are met:
#1. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#3. The name of the author may not be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
#IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
#IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
#NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""Local Simularity Analysis Package.

This python module provide tools for aligning and analyzing the time shift
dependent pairwise correlation between two sequence of evenly spaced 
observation data. Permutation test is used to estimate the P-value. 
The results can be summarized and easily queried for desired analysis.
"""
from setuptools import setup, find_packages
from distutils.core import Extension
from distutils.command import build
import os, sys

doclines=__doc__.splitlines()

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

#compcore_module = Extension('lsa/_compcore', sources=['lsa/compcore_wrap.cpp', 'lsa/compcore.cpp'] , )
# in above settings, put compcore.cpp, compcore.h, compcore.i file in lsa subdir
# use swig -python -c++ lsa/compcore.i -o compcore_wrap.cpp 
# then use python setup.py build
# _compcore.so shall appear in lsa subdir in build/build-arch dir
# when import module shall use import lsa.xxx

class my_build(build.build): 
    # different order: build_ext *before* build_py, so that 
    # build_py can already use ctypes! 
    sub_commands = [('build_ext', build.build.has_ext_modules), 
        ('build_py', build.build.has_pure_modules), 
        ('build_clib', build.build.has_c_libraries), 
        ('build_scripts', build.build.has_scripts), ] 

setup(name="lsa",
    version="1.0.0",
    description=doclines[0],
    long_description="\n".join(doclines[2:]),
    author="Li Charlie Xia",
    author_email="lxia@usc.edu",
    url="http://meta.usc.edu/softs/lsa",
    license="BSD",
    platforms=["Linux"],
    packages=find_packages(exclude=['ez_setup', 'test', 'doc']),
    include_package_data=True,
    #packages=['lsa'],
    #requires=['numpy(>=1.1)','scipy(>=0.6)','python(>=2.5)','matplotlib(>=0.98)'],
    zip_safe=False,
    install_requires=[ "python >= 2.7","numpy >= 1.0","scipy >= 0.6" ],
    provides=['lsa'],
    ext_modules = [Extension('lsa._compcore', sources = ['lsa/compcore_wrap.cpp', 'lsa/compcore.cpp'],
                                              depends = ['lsa/compcore.h'],
                   )],
    #ext_modules = [Extension('lsa._compcore', sources = ['lsa/compcore.i', 'lsa/compcore.cpp'],
    #                               depends=['lsa/compcore.h'],
    #                               swig_opts=['-c++', '-nomodern', '-classic', '-nomodernargs'])],
    py_modules = ['lsa.compcore', 'lsa.lsalib', 'lsa.lsaio'],
    #scripts = ['lsa-standalone/lsa-compute.py', 'lsa-standalone/lsa-query.py', 'lsa-standalone/lsa-infer.py'],
    cmdclass = {'build': my_build},
    data_files = [('',['INSTALL.txt','LICENSE.txt', ])],
    entry_points = { 
        'console_scripts': [
            'lsa_compute = lsa.lsa_compute:main',
            'lsa_query = lsa.lsa_query:main',
        ]
    }
)
