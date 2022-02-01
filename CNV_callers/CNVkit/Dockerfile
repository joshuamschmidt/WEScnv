FROM python:3.10

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

WORKDIR /usr/src/app

RUN pip install --no-cache-dir cnvkit

COPY . .

WORKDIR /data

CMD ["/bin/bash"]

