FROM python:3.7.8-alpine
LABEL maintainer="dirk.winkelhardt@gmail.com"

WORKDIR /usr/src/macpepdb

COPY . .
RUN apk update \
    && apk add --no-cache libxml2-dev libxslt-dev postgresql-dev gcc libc-dev\
    && pip install --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt \
    && pip install .



ENTRYPOINT [ "python", "-m", "macpepdb" ]
