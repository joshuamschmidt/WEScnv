FROM rocker/r-ver:4.1
# image for FNN
# rocker/r-ver is based on debian

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

# FNN and data.table
RUN install2.r -e  FNN
RUN install2.r -e  data.table
RUN install2.r -e stringdist
RUN install2.r -e cluster
RUN install2.r -e svglite
RUN install2.r -e ggplot2
RUN install2.r -e ggrepel
RUN install2.r -e factoextra
WORKDIR /data
