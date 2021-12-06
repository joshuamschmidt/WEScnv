FROM python:3

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

WORKDIR /tmp/install

COPY conifer_v0.2.2/ /tmp/install

CMD ["/bin/bash"]

RUN pip install tables

RUN chmod +x conifer.py && \
  mv conifer* /usr/local/bin && \
  cd && \
  rm -rf /tmp/install

WORKDIR /data




