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
import argparse
import datetime
import os
import pprint

script = '''
#!/bin/bash
#PBS -N geoglam_h{h}v{v}
#PBS -P {project} 
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l mem=64GB
#PBS -l ncpus=12
#PBS -o {log_dir}
#PBS -e {log_dir}

set -e
set -u
set -x

PYTHON={miniconda_path}/bin/python

cd {base_dir}
bash pipeline.sh {h} {v} {year_start} {year_end} fc_prod {prod_dir} $PYTHON 2>&1|tee {log_dir}/geoglam.h{h}v{v}.log 
bash pipeline.sh {h} {v} {year_start} {year_end} all_agg {prod_dir} $PYTHON 2>&1|tee -a {log_dir}/geoglam.h{h}v{v}.log 

'''

parser = argparse.ArgumentParser()
parser.add_argument('--project', type=str, default='', help='nci project')
parser.add_argument('--h', type=int, default=-1, help='h of tile')
parser.add_argument('--v', type=int, default=-1, help='v of tile')
parser.add_argument('--year_start', type=int, default=2001, help='staring year')
curr_year = datetime.datetime.now().year
parser.add_argument('--year_end', type=int, default=curr_year, help='ending year')
parser.add_argument('--prod_dir', type=str, default='', help='directory of the products')
parser.add_argument('--log_dir', type=str, default='', help='directory for pbs logs')
args = parser.parse_args()

pprint.pprint(args)

assert len(args.project) > 0, 'please specify NCI project'
assert args.h >= 0 and args.v >= 0
assert os.path.isdir(args.prod_dir), "product directory '%s' not found" % args.prod_dir
assert os.path.isdir(args.log_dir), "pbs log directory '%s' not found" % args.log_dir

curr_dir = os.path.dirname(os.path.realpath(__file__))
assert os.getcwd() == curr_dir, \
  "Please be sure that you have cd'ed into the directory of this python script"

miniconda_path = os.path.join(curr_dir, 'miniconda')
assert os.path.isdir(miniconda_path), \
  '''miniconda not found. Possible reasons: 
  1) Have you successfully run setup_env.sh?
  2) Is miniconda in the same directory of this python script?'''

h_str = '%02d' % args.h
v_str = '%02d' % args.v

pbs = script.format(h=h_str, v=v_str,
  miniconda_path=miniconda_path, year_start=args.year_start, base_dir=curr_dir,
  project=args.project, year_end=args.year_end, prod_dir=args.prod_dir, 
  log_dir=args.log_dir)

print pbs
with open('%s/h%sv%s.pbs.sh' % (args.log_dir, h_str, v_str), 'w') as fid:
  fid.write(pbs)

