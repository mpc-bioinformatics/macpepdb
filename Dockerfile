FROM python:3.9-buster
LABEL maintainer="dirk.winkelhardt@rub.de"

WORKDIR /usr/src/macpepdb

COPY . .
RUN apt-get update -y \
    && apt-get install -y libxml2-dev libxslt-dev libpq-dev gcc g++ libc-dev libev-dev \
    && pip install --upgrade pip \
    && pip install pipenv \
    && pipenv install --system


ENTRYPOINT [ "python", "-m", "macpepdb" ]
