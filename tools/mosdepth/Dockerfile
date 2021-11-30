# Dockerfile for mosdepth
FROM alpine:3.15

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apk update && \
    apk upgrade && \
    apk add bash

CMD ["/bin/bash"]

RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.2/mosdepth && \
  chmod 755 mosdepth && \
  mv mosdepth /usr/local/bin/

WORKDIR /data

