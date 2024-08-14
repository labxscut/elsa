FROM ubuntu:latest
# Set the working directory
WORKDIR /setup
# Install prerequisites
RUN apt-get update --fix-missing
RUN apt-get -y install curl git git-lfs build-essential                   
RUN apt-get -y install python3 python3-dev python3-full                               
RUN apt-get -y install python3-numpy python3-scipy python3-matplotlib 
RUN apt-get -y install python-is-python3 python-dev-is-python3 
RUN apt-get -y install pipx                              
# Install elsa
RUN git clone --verbose https://github.com/labxscut/elsa.git 
RUN cd elsa && git lfs checkout py3
RUN cd elsa && git lfs pull
RUN cd elsa && pipx install .
# Run elsa prompt
RUN lsa_compute --help
