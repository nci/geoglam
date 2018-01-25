# Copyright 2017 National Computational Infrastructure(NCI). 
# All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =========================================================================
#!/bin/bash
set -e
set -x 
set -u

NCI_PROJECT=<nci project number e.g. fr1>
BASE_DIR=<absoulte full path>/geoglam_data
DATA_DIR=$BASE_DIR/data
LOG_DIR=$BASE_DIR/logs

#cd monthly_medoids
#../miniconda/bin/python compile_medoids.py
#cd ..

mkdir -p $DATA_DIR
mkdir -p $LOG_DIR
job_hv_file=$LOG_DIR/job_hv.csv.`date +%Y-%m-%d_%H-%M-%S`
echo jobid,h,v > $job_hv_file

TILES=$(cat $1)
for tile in $TILES
do
  h=$(echo $tile|cut -d',' -f1)
  v=$(echo $tile|cut -d',' -f2)
  python pbs_geoglam.py --project $NCI_PROJECT --h $h --v $v --prod_dir $DATA_DIR --log_dir $LOG_DIR
  jobid=$(qsub $LOG_DIR/h${h}v${v}.pbs.sh)
  echo $jobid,$h,$v >> $job_hv_file
  rm $LOG_DIR/h${h}v${v}.pbs.sh

done

