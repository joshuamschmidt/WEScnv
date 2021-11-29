FROM rocker/r-ver:4.1

# rocker/r-ver is based on debian

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

# create a working directory and work from there
RUN mkdir /tmp/install
WORKDIR /tmp/install

RUN apt-get update && apt-get upgrade -y && \
	apt-get install -y curl gengetopt git libblas-dev libbz2-dev libcurl4-openssl-dev libc6-dev libncurses-dev lzma liblzma-dev perl python3-venv python3-pip wget zlib1g-dev 

RUN git clone https://github.com/ebiggers/libdeflate.git && \
	cd libdeflate && \
	make && \
	make install && \
	make clean

RUN wget -qO- https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2 | tar xvjf - && \
    cd htslib-1.14 && \
    make && \
    make install && \
    make clean

RUN wget -qO- https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 | tar xvjf - && \
    cd samtools-1.14 && \
    make && \
    make install && \
    make clean

RUN wget -qO- https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2 | tar xvjf - && \
    cd bcftools-1.14 && \
    make && \
    make install && \
    make clean

RUN wget -qO- https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v2.4.40.tar.bz2 | tar xvjf - && \
	cp bin/* /usr/local/bin

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
	mv bedtools.static.binary /usr/local/bin/bedtools && \
	chmod a+x /usr/local/bin/bedtools

RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.2/mosdepth && \
	chmod a+x mosdepth && \
	mv mosdepth /usr/local/bin/

# RUN wget https://github.com/brentp/hts-nim/releases/download/v0.2.8/hts_nim_static_builder && \
# 	chmod a+x hts_nim_static_builder && \
# 	mv hts_nim_static_builder /usr/local/bin/

RUN apt-get install -y liblapack-dev build-essential

RUN curl https://nim-lang.org/choosenim/init.sh -sSf > init.sh && \
	sh init.sh -y && \
	rm init.sh

ENV PATH=/root/.nimble/bin:$PATH
RUN git clone https://github.com/brentp/hts-nim.git && \
	cd hts-nim && \
	nimble install -y

RUN git clone https://github.com/brentp/hts-nim-tools.git && \
 	cd hts-nim-tools && \
 	nimble install -y && \
 	mv hts_nim_tools /usr/local/bin/

# CNV BINARIES
#XHMM
RUN git clone https://bitbucket.org/statgen/xhmm.git && \
	cd xhmm && \
	make CXXFLAGS='-Wno-error=deprecated-declarations' && \
	chmod a+x /tmp/install/xhmm/build/execs/xhmm && \
	mv /tmp/install/xhmm/build/execs/xhmm /usr/local/bin/

# CLAMMS
RUN git clone https://github.com/rgcgithub/clamms.git && \
	cd clamms && \
	make && \
	find . -maxdepth 1 -type f -executable | xargs -n1 cp -t /usr/local/bin

ENV CLAMMS_DIR=/usr/local/bin

# Tools in R
RUN install2.r -e data.table
RUN install2.r -e BiocManager
ENV BIOC /usr/local/lib/R/site-library/littler/examples

# ExomeDepth
RUN "$BIOC"/installBioc.r Biostrings IRanges Rsamtools GenomicRanges GenomicAlignments
RUN install2.r -e  ExomeDepth

# CODEX2
RUN apt-get install -y libxml2
RUN "$BIOC"/installBioc.r CODEX

# cn.MOPS
RUN "$BIOC"/installBioc.r cn.mops

RUN rm -rf ../downloaded_packages && \
	rm -rf /tmp

WORKDIR /app


# Tools in python
# python3
# CNVkit

# python 2.7 
#curl -sSL https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
# CONINFER
