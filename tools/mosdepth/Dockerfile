# Dockerfile for mosdepth
FROM alpine:3.15

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apk update && \
    apk upgrade && \
    apk add bash

RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.2/mosdepth && \
  chmod a+x mosdepth && \
  mv mosdepth /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/mosdepth"]
