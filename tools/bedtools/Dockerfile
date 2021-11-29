# Dockerfile for bedtools
FROM alpine:3.15

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apk update && \
  apk upgrade && \
  apk add bash

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
  mv bedtools.static.binary /usr/local/bin/bedtools && \
  chmod a+x /usr/local/bin/bedtools

WORKDIR /app


