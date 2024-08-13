FROM ubuntu:focal
# Set the working directory (will return to WORKDIR after each RUN)
WORKDIR /setup
# Install prerequisites
RUN apt-get update --fix-missing
RUN apt-get -y install curl git git-lfs build-essential                   # install curl, git and build tools
RUN apt-get -y install python2 python2-dev                             # install python2
RUN apt-get -y install python-is-python2 python-dev-is-python2 # set python2 to system python
RUN apt-get -y install python-setuptools                                  # install python2 and necessaries
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py | python2  # install pip
RUN pip install numpy scipy                                                    # install numpy and scipy
# Install elsa
RUN git lfs clone --verbose https://bitbucket.org/charade/elsa.git
RUN cd elsa && python2 setup.py install --force
# Run elsa
RUN lsa_compute --help
