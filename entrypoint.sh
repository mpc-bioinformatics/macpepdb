#!/bin/bash

. ~/appenv/bin/activate

# If macpepdb web serve is called with option gunicorn build params and pass it to gunicorn
# otherwise pass parameter to macpepdb.
if [[ "$1" == "web" ]] && [[ "$2" == "serve" ]] && [[ "$@" == *"--gunicorn"* ]]
then
    gunicorn_args="$(python -m macpepdb "$@")"
    set -- gunicorn "$gunicorn_args"
    eval "$@"
else
    set -- python -m macpepdb "$@"
    exec "$@"
fi