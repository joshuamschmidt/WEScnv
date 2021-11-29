# Dockerfile for bedops
FROM alpine:3.15

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apk update && \
  apk upgrade && \
  apk add bash

RUN wget -qO- https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v2.4.40.tar.bz2 | tar xvjf - && \
  cp bin/* /usr/local/bin

WORKDIR /app
