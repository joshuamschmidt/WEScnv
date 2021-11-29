# Dockerfile for CLAMMS
FROM alpine:3.15 AS builder

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apk update && \
  apk add bash && \
  apk add gcc && \
  apk add g++ && \
  apk add make && \
  apk add curl && \
  apk add zlib && \
  apk add zlib-dev && \
  apk add libbz2 && \
  apk add xz-dev && \
  apk add libcurl && \
  apk add libressl && \
  apk add lapack && \
  apk add lapack-dev && \
  apk add libpthread-stubs && \
  apk add gengetopt && \
  apk add git && \
  apk add nim && \
  apk add perl && \
  apk add python3 && \
  apk add R && \
  apk add bzip2-dev && \
  apk add curl-dev && \
  apk add ncurses-dev

WORKDIR /tmp

RUN git clone https://github.com/rgcgithub/clamms.git && \
  cd clamms && \
  make && \
  find . -maxdepth 1 -type f -executable | xargs -n1 cp -t /usr/local/bin

FROM alpine:3.15

RUN apk update && \
  apk upgrade && \
  apk add bash

COPY --from=builder /usr/local/bin /usr/local/bin/
WORKDIR /app
ENV CLAMMS_DIR=/usr/local/bin
