#NOTE: developer tested the commands by folloing actions:
#      docker pull ubuntu:latest
#      docker images
#      docker run --memory=2g -i -t ubuntu:latest /bin/bash
#      docker run --memory=2g -i -t id /bin/bash
#      docker-machine scp Dockerfile main:~/zoomx
#      docker build --memory=2g ~/zoomx
#      docker start --memory=2g 301176b69086
#      docker exec -it 301176b69086 /bin/bash
#NOTE: merging RUN as file become stable as every RUN creates a commit which has a limit

##### Following Replaced By from charade/xlibbox:basic
#FROM ubuntu:xenial
#
#MAINTAINER Charlie Xia <xia.stanford@gmail.com>
#
#WORKDIR $HOME/
#RUN ulimit -s unlimited
#RUN mkdir $HOME/setup && mkdir $HOME/bin
#
#### install git and dev files ###
##NOTE: RUN echo "PATH=$PATH:$HOME/bin" >>/etc/enviroment #not working, why?
#RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
#RUN gpg -a --export E084DAB9 >$HOME/rstudio.key
#RUN apt-key add $HOME/rstudio.key
#RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" >>/etc/apt/sources.list
#RUN echo "PATH=$PATH:$HOME/bin" >>/etc/enviroment
#RUN apt-get update
#RUN apt-get install -y build-essential software-properties-common
#RUN apt-get install -y bzip2 libbz2-dev liblzma-dev openssl libssl-dev
#RUN apt-get install -y zlib1g-dev libncurses5-dev wget git unzip libcurl4-openssl-dev libxml2-dev
#
#### python, has to precede bedtools ###
#RUN apt-get install -y python libpython-dev python-pip
#RUN apt-get install -y python-numpy python-scipy
#### install python packages ###
#RUN pip install -U numpy scipy tables six pandas pysam pybedtools dendropy
#
#### install bedtools ###
#RUN cd $HOME/setup && git clone https://github.com/arq5x/bedtools.git
#RUN cd $HOME/setup/bedtools && make && cp bin/* $HOME/bin
#
#### install samtools ###
#RUN cd $HOME/setup && git clone https://github.com/samtools/samtools.git
#RUN cd $HOME/setup && git clone https://github.com/samtools/htslib.git
#RUN cd $HOME/setup/samtools && make && cp samtools $HOME/bin
##### Above Replaced By from charade/xlibbox:basic

FROM charade/xlibbox:basic
RUN apt-get update

### avoid RPC error with https ###
RUN git config --global http.sslVerify false
RUN git config --global http.postBuffer 1048576000

### install elsa ###
RUN cd $HOME/setup && git clone --verbose https://charade@bitbucket.org/charade/elsa.git
RUN cd $HOME/setup/elsa && python setup.py install --force

### run test scripts ###
RUN cd $HOME/setup/elsa/test && ./test.sh
