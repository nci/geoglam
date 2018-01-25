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
import re
import glob
import argparse

parser = argparse.ArgumentParser(description="""Parse Geoglam logs for errors""")
parser.add_argument("log_dir", default="", help="Director for log files")    
parser.add_argument("--output_file", default="failed_tiles.csv", help="Output file")    
args = parser.parse_args()

regex = r"\+ echo __checkpoint_cleanup"
pattern = re.compile(regex)

log_files = glob.glob(args.log_dir + '/geoglam.*.log')
print '%d log files found' % len(log_files)

failed_tiles = ''
for i, f in enumerate(log_files):
  with open(f, 'r') as fid:
    log = fid.read()
  matches = re.findall(regex, log)

  msg = '(%d/%d)' % (i+1, len(log_files))
  if len(matches) != 1:
    filename = f.split('/')[-1]
    hv = filename.split('.')[1]
    parts = hv.split('v')
    v = parts[1]
    h = parts[0][1:]
    failed_tiles += '%s,%s' % (h, v)
    
    if i < len(log_files) - 1:
      failed_tiles += '\n'

    msg += 'FAIL %s' % f
  else:
    msg += 'OK %s' % f

  print msg
    
if len(failed_tiles) > 0:
  print 'Saving failed tiles to %s' % args.output_file
  with open(args.output_file, 'w') as fid:
    fid.write(failed_tiles)
else:
  print 'No failed tiles found'
  

