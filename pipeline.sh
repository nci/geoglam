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

date "+%Y-%m-%d %H:%M:%S"
VER=310

H=$1
V=$2
YEAR_BGN=$3
YEAR_END=$4
CHECKPOINT=$5

BASE_OUTPUT=$6/v${VER}/tiles
PYTHON=$7

PARALLEL=utils/parallel
FC_OUTPUT=$BASE_OUTPUT/8-day/cover/fc_outputs_tmp
FC_LOGS=$BASE_OUTPUT/8-day/cover/fc_logs
PACKED_OUTPUT=$BASE_OUTPUT/8-day/cover
MEDOIDS_OUTPUT=$BASE_OUTPUT/monthly/cover
ANOMALY_MEAN_DIFF_OUTPUT=$BASE_OUTPUT/monthly/anomalies
ANOMALY_MEAN_DIFF_TMP_OUTPUT=${ANOMALY_MEAN_DIFF_OUTPUT}_h${H}v${V}_tmp
ANOMALY_PERCENTILE_OUTPUT=$BASE_OUTPUT/monthly/deciles
ANOMALY_PERCENTILE_TMP_OUTPUT=${ANOMALY_PERCENTILE_OUTPUT}_h${H}v${V}_tmp

if [[ "$CHECKPOINT" = 'fc_prod' ]]; then
mkdir -p $FC_OUTPUT
mkdir -p $PACKED_OUTPUT
mkdir -p $FC_LOGS

for i in $(seq $YEAR_BGN $YEAR_END)
do
rm -rf $FC_OUTPUT/h${H}v${V}.$i
mkdir -p $FC_OUTPUT/h${H}v${V}.$i

FREEMEM=$(free -k|grep 'cache:'|awk '{print $4}')
N_JOBS=$(echo "$FREEMEM/(4.5*1024^2)"|bc)

UPDATE_FILE=$(python fc_prod/get_update_file.py $PACKED_OUTPUT $i $H $V --version $VER)
seq 1 12|$PARALLEL -j$N_JOBS -k --joblog $FC_LOGS/h${H}v${V}.$i.parjob.log \
     $PYTHON fc_prod/main.py /g/data2/u39/public/data/modis/lpdaac-tiles-c6/ $i {} $H $V $FC_OUTPUT \
     --update_file $UPDATE_FILE --version $VER
$PYTHON fc_prod/packer.py $FC_OUTPUT/h${H}v${V}.$i/ $PACKED_OUTPUT --update_file $UPDATE_FILE --version $VER
rm -rf $FC_OUTPUT/h${H}v${V}.$i

date "+%Y-%m-%d %H:%M:%S"
done

date "+%Y-%m-%d %H:%M:%S"
echo __checkpoint_fc_prod
fi


if [[ "$CHECKPOINT" = 'all_agg' ]] || [[ "$CHECKPOINT" = 'medoids' ]]; then
FREEMEM=$(free -k|grep 'cache:'|awk '{print $4}')
N_JOBS=$(echo "$FREEMEM/(14.5*1024^2)"|bc)
seq $YEAR_BGN $YEAR_END|$PARALLEL -j$N_JOBS -k $PYTHON monthly_medoids/main.py $PACKED_OUTPUT $H $V {} $MEDOIDS_OUTPUT --version $VER

date "+%Y-%m-%d %H:%M:%S"
echo __checkpoint_medoids
fi


INIT_YEAR=2001
CURR_YEAR=$(date "+%Y")
if [[ "$CHECKPOINT" = 'all_agg' ]] || [[ "$CHECKPOINT" = 'anomaly_detection' ]]; then
FREEMEM=$(free -k|grep 'cache:'|awk '{print $4}')
N_JOBS=$(echo "$FREEMEM/(7.5*1024^2)"|bc)
mkdir -p ${ANOMALY_MEAN_DIFF_TMP_OUTPUT} 
mkdir -p ${ANOMALY_PERCENTILE_TMP_OUTPUT} 
seq 1 12|$PARALLEL -j$N_JOBS -k $PYTHON anomaly_detection/main.py $PACKED_OUTPUT $MEDOIDS_OUTPUT $H $V \
    $INIT_YEAR $CURR_YEAR {} ${ANOMALY_MEAN_DIFF_TMP_OUTPUT} ${ANOMALY_PERCENTILE_TMP_OUTPUT} --version $VER

date "+%Y-%m-%d %H:%M:%S"
echo __checkpoint_anomaly_detection
fi


if [[ "$CHECKPOINT" = 'all_agg' ]] || [[ "$CHECKPOINT" = 'combine_mean_diff' ]]; then
echo 'combine monthly mean differences'
$PYTHON anomaly_detection/combine_outputs.py $PACKED_OUTPUT ${ANOMALY_MEAN_DIFF_TMP_OUTPUT} $H $V \
    $INIT_YEAR $CURR_YEAR $ANOMALY_MEAN_DIFF_OUTPUT Mean_Diff f4 --version $VER

date "+%Y-%m-%d %H:%M:%S"
echo __checkpoint_combine_mean_diff
fi


if [[ "$CHECKPOINT" = 'all_agg' ]] || [[ "$CHECKPOINT" = 'combine_percentiles' ]]; then
echo 'combine monthly percentiles'
$PYTHON anomaly_detection/combine_outputs.py $PACKED_OUTPUT ${ANOMALY_PERCENTILE_TMP_OUTPUT} $H $V \
    $INIT_YEAR $CURR_YEAR $ANOMALY_PERCENTILE_OUTPUT Percentile u1 --version $VER

date "+%Y-%m-%d %H:%M:%S"
echo __checkpoint_combine_percentiles
fi


if [[ "$CHECKPOINT" = 'all_agg' ]] || [[ "$CHECKPOINT" = 'cleanup' ]]; then
rm -rf ${ANOMALY_MEAN_DIFF_TMP_OUTPUT}
rm -rf ${ANOMALY_PERCENTILE_TMP_OUTPUT}

date "+%Y-%m-%d %H:%M:%S"
echo __checkpoint_cleanup
fi

date "+%Y-%m-%d %H:%M:%S"




