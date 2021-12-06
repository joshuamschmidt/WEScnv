FROM rocker/r-ver:4.1
# image for CODEX2 1.26.0
# rocker/r-ver is based on debian

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apt-get update && apt-get upgrade -y && \
  apt-get install -y libbz2-dev libcurl4-openssl-dev liblzma-dev libxml2 zlib1g-dev

RUN install2.r -e data.table stringdist

RUN install2.r -e BiocManager

ENV BIOC /usr/local/lib/R/site-library/littler/examples

# ExomeDepth
RUN "$BIOC"/installBioc.r CODEX

CMD ["/bin/bash"]

COPY CODEX2.R /usr/local/bin

RUN chmod +x /usr/local/bin/CODEX2.R
