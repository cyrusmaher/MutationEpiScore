# EpiScore
Scoring for predicting spread of SARS-CoV-2 mutations

# Installation
## Using conda
The python environment can be configured with:
`conda env create -f environment.yml`

## Using docker
You can build the docker environment with
`docker build .`

An example for running from docker can be found in `run_docker.sh`
You will have to update `CONTAINER_ID` on the first line of the script to match the ID of the container you built.

# Usage
`python forecasting.py [GISAID metadata file] [Output folder]`

# Output
A table of EpiScores and EpiScore components for each observed mutation


