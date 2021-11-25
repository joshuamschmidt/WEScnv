FROM alpine:3.15
MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apk update && \
	apk add bash && \
	apk add gcc && \
	apk add g++ && \
	apk add make && \
	apk add curl && \
	apk add zlib && \
	apk add zlib-dev && \
	apk add libbz2 && \
	apk add xz-dev && \
	apk add libcurl && \
	apk add libressl && \
	apk add lapack && \
	apk add libpthread-stubs && \
	apk add gengetopt && \
	apk add git && \
	apk add nim && \
	apk add perl && \
	apk add python3 && \
	apk add R && \
	apk add bzip2-dev && \
	apk add curl-dev && \
	apk add ncurses-dev

WORKDIR /tmp

RUN wget -qO- https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2 | tar xvjf - && \
    cd htslib-1.14 && \
    make && \
    make install

RUN wget -qO- https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 | tar xvjf - && \
    cd samtools-1.14 && \
    make && \
    make install

RUN wget -qO- https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2 | tar xvjf - && \
    cd bcftools-1.14 && \
    make && \
    make install

RUN wget -qO- https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v2.4.40.tar.bz2 | tar xvjf - && \
	cp bin/* /usr/local/bin

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
	mv bedtools.static.binary /usr/local/bin/bedtools && \
	chmod a+x /usr/local/bin/bedtools

RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.2/mosdepth && \
	chmod a+x mosdepth && \
	mv mosdepth /usr/local/bin/

RUN wget https://github.com/brentp/hts-nim/releases/download/v0.2.8/hts_nim_static_builder && \
	chmod a+x hts_nim_static_builder && \
	mv hts_nim_static_builder /usr/local/bin/

RUN apk add nimble && \
	apk add lapack-dev

RUN git clone https://github.com/brentp/hts-nim.git && \
	cd hts-nim && \
	nimble install -y

RUN git clone https://github.com/brentp/hts-nim-tools.git && \
 	cd hts-nim-tools && \
 	nimble install -y && \
 	mv hts_nim_tools /usr/local/bin/

# CNV BINARIES
# XHMM
RUN git clone https://bitbucket.org/statgen/xhmm.git && \
	cd xhmm && \
	make CXXFLAGS='-Wno-error=deprecated-declarations' && \
	chmod a+x /tmp/xhmm/build/execs/xhmm && \
	mv /tmp/xhmm/build/execs/xhmm /usr/local/bin/

# CLAMMS
RUN git clone https://github.com/rgcgithub/clamms.git && \
	cd clamms && \
	make && \
	find . -maxdepth 1 -type f -executable | xargs -n1 cp -t /usr/local/bin

ENV CLAMMS_DIR=/usr/local/bin

# CNVkit

# ExomeDepth
RUN R -e "install.packages('ExomeDepth',dependencies=TRUE, repos='http://cran.rstudio.com/',lib='/usr/lib/R/library')"

# Wise