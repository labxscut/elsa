# From Ubuntu Noble
FROM ubuntu:latest
# All Necessary configs
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
WORKDIR /setup

# Install C prerequisites
RUN apt-get update --fix-missing && \
    apt-get -y install curl git git-lfs build-essential && \
    rm -rf /var/lib/apt/lists/*

# Install Python and Conda
RUN apt-get update --fix-missing && \
    apt-get -y install python3 python3-dev python3-pip && \
    curl -sLo miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash miniconda.sh -b -p ~/miniconda && \
    echo 'export PATH=~/miniconda/bin:$PATH' >> ~/.bashrc && \
    rm miniconda.sh && \
    rm -rf /var/lib/apt/lists/*

# Install elsa
ENV PATH=~/miniconda/bin:$PATH
RUN conda create -n elsa -y python=3.8 && \ 
    echo "source activate elsa" >> ~/.bashrc && \
    . ~/miniconda/etc/profile.d/conda.sh && \
    conda activate elsa && \
    pip install numpy scipy matplotlib && \
    git clone https://github.com/labxscut/elsa.git && \
    cd elsa && \
    git checkout devel && \
    git pull && \
    pip install . && \
    lsa_compute --help && \
    cd test && \
    . test.sh
