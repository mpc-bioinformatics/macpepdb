FROM mambaorg/micromamba:1.2.0-jammy
LABEL maintainer="dirk.winkelhardt@rub.de"

ARG NEW_MAMBA_USER_ID=1000
ARG NEW_MAMBA_USER_GID=1000

USER root
RUN usermod "--login=mambauser" "--home=/home/mambauser" \
        --move-home "-u ${NEW_MAMBA_USER_ID}" "${MAMBA_USER}" && \
    groupmod "--new-name=mambauser" \
        "-g ${NEW_MAMBA_USER_GID}" "${MAMBA_USER}" && \
    # Update the expected value of MAMBA_USER for the
    # _entrypoint.sh consistency check.
    echo "mambauser" > "/etc/arg_mamba_user" && \
    :
RUN apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y git curl

WORKDIR /home/mambauser
COPY --chown=mambauser:mambauser . macpepdb/
WORKDIR /home/mambauser/macpepdb

USER mambauser
ENV HOME /home/mambauser
ENV PATH $PATH:$HOME/.local/bin
ENV ENV_NAME=macpepdb

RUN echo 'show_banner: false' > ~/.mambarc

RUN micromamba env create -y -f environment.yml \
    && micromamba clean --all --yes


ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "/home/mambauser/macpepdb/entrypoint.sh" ]
