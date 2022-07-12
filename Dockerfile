FROM python:3.9-buster as app
LABEL maintainer="dirk.winkelhardt@rub.de"

ARG USER_ID=999
ARG GROUP_ID=999

ENV USER_ID=$USER_ID
ENV GROUP_ID=$GROUP_ID

WORKDIR /home/app
# Copy macpepdb
COPY . macpepdb_src/

RUN apt-get update -y && apt-get install -y libxml2-dev libxslt-dev libpq-dev gcc g++ libc-dev \
    && groupadd -g $GROUP_ID app \
    && useradd -g $GROUP_ID -m -s /bin/bash -u $USER_ID app \
    && chown -R app:app .

USER app
ENV HOME /home/app
ENV PATH $PATH:$HOME/.local/bin

RUN python -m pip install --user virtualenv \
    && virtualenv ~/appenv \
    && . ~/appenv/bin/activate \
    && pip install --upgrade pip \
    && pip install ~/macpepdb_src

ENTRYPOINT [ "/home/app/macpepdb_src/entrypoint.sh" ]
