FROM ubuntu:latest
# Set the working directory
WORKDIR /setup
# Install prerequisites
RUN apt-get update --fix-missing
RUN apt-get -y install curl git git-lfs build-essential                   
RUN apt-get -y install python3 python3-dev                                
RUN apt-get -y install python3-numpy python3-scipy python3-matplotlib 
RUN apt-get -y install python-is-python3 python-dev-is-python3 
RUN apt-get -y install pip                                
# Install elsa
RUN git lfs clone --verbose https://github.com/labxscut/elsa.git
RUN cd elsa && python setup.py install --force
# Run elsa prompt
RUN lsa_compute --help
