FROM rocker/r-ver:4.1
# image for ExomeDepth 1.1.15
# rocker/r-ver is based on debian

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN install2.r -e data.table pacman
RUN install2.r -e BiocManager
ENV BIOC /usr/local/lib/R/site-library/littler/examples

RUN apt-get update && apt-get upgrade -y && \
  apt-get install -y libbz2-dev libcurl4-openssl-dev liblzma-dev  zlib1g-dev
# ExomeDepth
RUN "$BIOC"/installBioc.r Biostrings IRanges Rsamtools GenomicRanges GenomicAlignments
RUN install2.r -e  ExomeDepth

WORKDIR /app
