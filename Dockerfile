FROM python:3.8.8-alpine
LABEL maintainer="dirk.winkelhardt@gmail.com"

WORKDIR /usr/src/macpepdb

COPY . .
RUN apk update \
    && apk add --no-cache libxml2-dev libxslt-dev postgresql-dev gcc g++ libc-dev\
    && pip install --upgrade pip \
    && pip install pipenv \
    && pipenv install --system



ENTRYPOINT [ "python", "-m", "macpepdb" ]
