
CONTAINER_ID=0b8f16880f97  # The container ID returned after running `docker build .`
CODE=$(pwd) # Assuming you're running this script from the current directory
DATA=${HOME}/Data # Update to the directory where you want to store the data
OUT_DIR="forecast" # the name of the output directory within $DATA
METADATA_FNAME="metadata_oct2021.tsv"  # name of your GISAID metadata file

docker run -it \
-v ${CODE}:/root/Code \
-v ${DATA}:/root/Data \
${CONTAINER_ID} \
bash -c "conda run -n sars2-forecasting python /root/Code/forecasting.py  \
/root/Data/${METADATA_FNAME} \
/root/Data/${OUT_DIR}"
