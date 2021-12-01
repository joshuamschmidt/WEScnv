FROM python:3.10

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

WORKDIR /data

CMD ["/bin/bash"]

COPY countsToMatrix.py /usr/local/bin

RUN chmod +x /usr/local/bin/countsToMatrix.py
