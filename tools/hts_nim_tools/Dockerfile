# Dockerfile for hts_nim_tools
FROM alpine:3.15

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apk update && \
 apk upgrade && \
 apk add bash

COPY hts_nim_tools /usr/local/bin/
WORKDIR /data

