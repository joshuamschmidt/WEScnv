FROM python:3.10

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

WORKDIR /usr/src/app

RUN pip install --no-cache-dir cnvkit

COPY . .

WORKDIR /data

CMD ["/bin/bash"]

ENV MPLCONFIGDIR=/data

FROM rocker/r-ver:4.1
RUN apt-get update && apt-get upgrade -y && \
  apt-get install -y libbz2-dev libcurl4-openssl-dev liblzma-dev  zlib1g-dev

RUN install2.r -e BiocManager
ENV BIOC /usr/local/lib/R/site-library/littler/examples
RUN "$BIOC"/installBioc.r DNAcopy
