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
set -u
set -x

wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh

INSTALL_DIR=$(pwd)
bash miniconda.sh -b -p $INSTALL_DIR/miniconda 
rm miniconda.sh

export PATH=$INSTALL_DIR/miniconda/bin:$PATH 
which python
which conda

conda install -y numpy
conda install -y scipy
conda install -y cython 
conda install -y -c anaconda pillow

conda install -y -c anaconda netcdf4 
conda install -y -c anaconda gdal=2.2.2

mkdir -p ~/.parallel
touch ~/.parallel/will-cite
echo 'testing GNU parallel...'
seq 1 2|utils/parallel echo {}

echo 'testing environment setup...'
python test_env.py

echo 'compiling fast_medoids.pyx...'
cd monthly_medoids
python compile_medoids.py

echo 'Environment successfully set up'



