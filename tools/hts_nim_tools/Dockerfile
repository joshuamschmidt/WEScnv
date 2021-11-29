# Dockerfile for hts_nim_tools
FROM alpine:3.15 AS builder

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN apk update && \
  apk upgrade && \
  apk add bash && \
  apk add bzip2-dev && \
  apk add curl && \
  apk add curl-dev && \
  apk add g++ && \
  apk add gcc && \
  apk add git && \
  apk add lapack && \
  apk add lapack-dev && \
  apk add libbz2 && \
  apk add libcurl && \
  apk add libressl && \
  apk add make && \
  apk add ncurses-dev && \
  apk add nim && \
  apk add nimble && \
  apk add perl && \
  apk add xz-dev && \
  apk add zlib && \
  apk add zlib-dev

WORKDIR /tmp/install

RUN wget -qO- https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2 | tar xvjf - && \
    cd htslib-1.14 && \
    make && \
    make install

RUN git clone https://github.com/brentp/hts-nim.git && \
  cd hts-nim && \
  nimble install -y

RUN git clone https://github.com/brentp/hts-nim-tools.git && \
  cd hts-nim-tools && \
  nimble install -y && \
  mv hts_nim_tools /usr/local/bin/

FROM alpine:3.15

RUN apk update && \
  apk upgrade && \
  apk add bash

COPY --from=builder /usr/local/bin/* /usr/local/bin/
COPY --from=builder /lib/ /lib/
COPY --from=builder /usr/local/lib /usr/local/lib/
COPY --from=builder /usr/local/include/htslib /usr/local/include/

WORKDIR /app

