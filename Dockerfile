FROM ubuntu:14.04
MAINTAINER Tom Whalley <twhalley93@gmail.com>

ENV PYMOL_VERSION 1.8.2.0
ENV HMMER_VERSION 3.2.1

#create bash
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

#basics
RUN apt-get update 
RUN apt-get install -y sudo wget curl make python python-pip pkg-config

#pymol specific dependencies
RUN apt-get install -y freeglut3 freeglut3-dev libpng3 libpng-dev \
libfreetype6 libfreetype6-dev pmw python-dev glew-utils libglew-dev \
libxml2-dev    libatlas-base-dev libgsl0-dev libblas-dev liblapack-dev \
gfortran libzmq1 libzmq-dev libc-dev   libtiff4-dev libjpeg8-dev \
zlib1g-dev liblcms2-dev libwebp-dev tcl8.5-dev tk8.5-dev python-tk

#python libraries
RUN pip install --upgrade pip
RUN pip install --ignore-installed six
RUN pip install numpy
RUN pip install scipy
RUN pip install biopython
RUN pip install matplotlib
RUN pip install argparse
RUN pip install swalign
RUN pip install networkx
RUN pip install pandas

#clean temp
RUN apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#install pymol
RUN wget --no-verbose https://sourceforge.net/projects/pymol/files/pymol/1.8/pymol-v${PYMOL_VERSION}.tar.bz2
RUN tar jxf pymol-v${PYMOL_VERSION}.tar.bz2
RUN rm pymol-v*
WORKDIR pymol
RUN python setup.py build install
WORKDIR /

#install HMMER
RUN wget http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz
RUN tar xvf hmmer-${HMMER_VERSION}.tar.gz
WORKDIR hmmer-${HMMER_VERSION}
RUN ./configure --enable-portable-binary --prefix=/hmmer && make && make install
RUN cd easel
RUN make install
ENV PATH="/hmmer/bin:${PATH}"
WORKDIR /

#ANARCI
COPY tools/anarci-1.3.tar.gz /
RUN tar -xvf anarci-1.3.tar.gz
WORKDIR anarci-1.3
RUN python setup.py install

#now CCP4
WORKDIR /
COPY tools/ccp4-7.0 /ccp4-7.0

#add the code
copy STACEI.py /
copy bin bin/

#activate shell
RUN echo "Activating CCP4 shell"
RUN source /ccp4-7.0/bin/ccp4.setup-sh

RUN echo "source ccp4-7.0/bin/ccp4.setup-sh" >> /root/.bashrc
RUN source /root/.bashrc