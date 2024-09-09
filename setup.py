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
import os, sys, subprocess

doclines=__doc__.splitlines()

#os.environ['CC'] = 'g++'  #temporary measure to trick distutils use g++, need update to distutils2
#lines = open("VERSION.txt", 'r').readlines()
#version_desc = ','.join([lines[1].strip(), lines[0].strip()])
print("works with python setup.py install or pipx install .", file=sys.stderr)

print("testing git availability ...", file=sys.stderr)
git_on_cmd = "echo 'def main():\n\t print('\\\"'$(cat VERSION.txt); @GIT: $(git log --pretty=format:'%h' | head -n 1)')\\\" > lsa/lsa_version.py"
try:
    subprocess.check_call(git_on_cmd, shell=True)
    print("Git commit number included in version info.", file=sys.stderr)
except subprocess.CalledProcessError:
    print("Git not available. Skipping commit number in version info.", file=sys.stderr)
    with open('lsa/lsa_version.py', 'w') as f:
        f.write("def main():\n\tprint('{}')".format(open('VERSION.txt').read().strip()))

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

#compcore_module = Extension('lsa/_compcore', sources=['lsa/compcore_wrap.cpp', 'lsa/compcore.cpp'] , )
# in above settings, put compcore.cpp, compcore.h, compcore.i file in lsa subdir
# and generate compcore.py file
# swig -python -c++ -o compcore_wrap.cpp compcore.i
# then use python setup.py build
# _compcore.so shall appear in lsa subdir in build/build-arch dir
# when import module shall use import lsa.xxx

class my_build(build.build):
    """
    Custom build class that changes the order of build sub-commands.

    This class modifies the order of build sub-commands to ensure that
    build_ext is executed before build_py. This allows build_py to use
    ctypes that may be generated during the build_ext phase.

    Attributes:
        sub_commands (list): A list of tuples representing the build
                             sub-commands in the desired order.
    """

    sub_commands = [
        ('build_ext', build.build.has_ext_modules),
        ('build_py', build.build.has_pure_modules),
        ('build_clib', build.build.has_c_libraries),
        ('build_scripts', build.build.has_scripts),
    ]

setup(name="lsa",
    version="2.0.0",
    description=doclines[0],
    long_description="\n".join(doclines[2:]),
    author="Li Charlie Xia",
    author_email="lcxia@scut.edu.cn",
    url="http://github.com/labxscut/elsa",
    license="BSD",
    platforms=["Linux"],
    packages=find_packages(exclude=['ez_setup', 'test', 'doc']),
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'numpy>=1.20.0',  # Updated to a version that supports the new float dtypes
        'scipy>=1.6.0',  # Ensure this version or higher is available
        'matplotlib>=3.3.0',
        # Add other dependencies as needed
    ],
    provides=['lsa'],
    ext_modules = [Extension('lsa._compcore', sources = ['lsa/compcore_wrap.cpp', 'lsa/compcore.cpp'],
                                              depends = ['lsa/compcore.hpp'],
                   )],
    #ext_modules = [Extension('lsa._compcore', sources = ['lsa/compcore.i', 'lsa/compcore.cpp'],
    #                               depends=['lsa/compcore.hpp'],
    #                               swig_opts=['-c++', '-nomodern', '-classic', '-nomodernargs'])],
    py_modules = ['lsa.compcore', 'lsa.lsalib_core', 'lsa.lsalib_stats', 'lsa.lsalib_normalization', 'lsa.lsalib_analysis', 'lsa.lsalib_utils', 'lsa.lsaio'],
    cmdclass = {'build': my_build},
    data_files = [('',['README.rst','LICENSE.txt','VERSION.txt'])],
    entry_points = { 
        'console_scripts': [
            'lsa_compute = lsa.lsa_compute:main',
            'lsa_query = lsa.lsa_query:main',
            'lsa_infer = lsa.lsa_infer:main',
            'lsa_sim = lsa.lsa_sim:main',
            'lsa_totrend = lsa.lsa_totrend:main',
            'lsa_para = lsa.lsa_para:main',
			'lsa_chkdat = lsa.lsa_chkdat:main',
			'lsa_fixqv = lsa.lsa_fixqv:main',
			'lsa_version = lsa.lsa_version:main'
        ]
    },
)
