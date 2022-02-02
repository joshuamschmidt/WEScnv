FROM ubuntu:latest

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base python-is-python3 python3-pip python3-setuptools python3-dev

WORKDIR /app

COPY requirements.txt /app/requirements.txt

RUN pip3 install -r requirements.txt

COPY . /app

RUN Rscript -e "install.packages('BiocManager')"

ENV BIOC /usr/local/lib/R/site-library/littler/examples

RUN Rscript -e "BiocManager::install('DNAcopy')"

ENV MPLCONFIGDIR=/app

COPY . /app
