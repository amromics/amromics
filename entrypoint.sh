#!/usr/bin/env bash

set -e

echo " "
echo "############### Docker Entrypoint: amromics Startup ###############"
echo " "
date

echo "##### Variables #####"
echo "PATH: $PATH"

echo "##### activate the micromamba env #####"
eval "$(micromamba shell hook --shell bash)" && 
micromamba activate amromics


echo " "
echo "############### Docker CMD ###############"
echo " "
exec "$@"
