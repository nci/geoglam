#!/bin/bash
set -xeu

METADATA_FILE="$1"
UPDATE_DIR="/g/data2/tc43/modis-fc/v310/tiles/"
PYTHON="../miniconda/bin/python"
PARALLEL="../utils/parallel"

LOG_DIR=logs
mkdir -p $LOG_DIR
LOG_FILE=$LOG_DIR/`date "+%Y-%m-%d_%H-%M-%S"`.log

time find "$UPDATE_DIR" -type f -name "*.nc" | $PARALLEL --will-cite --progress -j4 -k --joblog $LOG_FILE \
     $PYTHON update_metadata.py --update_file {} --metadata_file $METADATA_FILE
