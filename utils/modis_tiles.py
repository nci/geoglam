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
import os
import glob

base_dir = '/g/data2/u39/public/data/modis/lpdaac-tiles-c6/MCD43A4.006/'

dirs = os.listdir(base_dir)
hv_lookup = {}
for i, d in enumerate(dirs):
  path = os.path.join(base_dir, d)
  files = glob.glob(os.path.join(path, '*.hdf'))
  for f in files:
    if not os.path.isfile(f):
      continue
    filename = f.split('/')[-1]
    parts = filename.split('.')
    hv = parts[-4]
    assert len(hv) == 6

    h = hv[1:3]
    v = hv[4:]
    
    assert int(h) >= 0 and int(h) <= 36
    assert int(v) >= 0 and int(v) <= 18

    if hv not in hv_lookup:
      hv_lookup[hv] = [h, v, 1]    
    else:
      hv_lookup[hv][-1] += 1 

  if i % (len(dirs) / 10) == 0:
    print '%d of %d done' % (i, len(dirs))

with open('modis_tiles.csv', 'w') as fid:
  for k, v in hv_lookup.iteritems():
    fid.write('%s,%s,%d\n' % (v[0], v[1], v[2]))





