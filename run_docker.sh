
CONTAINER_ID=0b8f16880f97  # The container ID returned after running `docker build .`
CODE=$(pwd)
DATA=${HOME}/Data # Update to the directory where you want to store the data
METADATA_FNAME="metadata_oct2021.tsv"  # name of your GISAID metadata file

docker run -it -v ${CODE}:/root/Code \
-v ${DATA}:/root/Data \
${CONTAINER_ID} \
bash -c "conda run -n sars2-forecasting python ~/Code/forecasting.py  \
${DATA}/${METADATA_FNAME}  \
${DATA}/forecast"
