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
You will have to update `CONTAINER_ID` on the first line of the script to match the ID of the container you built. The comments within the script will make it clear where to update the input file, output directory, etc.

# Obtaining input data
The input to this method is GISAID metadata. You can apply for access to GISAID here:
https://www.gisaid.org/registration/register/

# Usage
`python forecasting.py [GISAID metadata file] [Output folder]`

# Output
A table of EpiScores and EpiScore components for each observed mutation


