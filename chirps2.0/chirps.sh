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

if [ "$#" -ne 1 ]
then
  echo "A year has to be provided. Usage: chirps YYYY"
  exit 1
fi

module purge
module load gdal/2.0.0
module load netcdf/4.3.3.1

rm chirps-v2.0.$1.dekads.nc
wget ftp://ftp.chg.ucsb.edu/pub/org/chg/products/CHIRPS-2.0/global_dekad/netcdf/chirps-v2.0.$1.dekads.nc
nccopy -3 chirps-v2.0.$1.dekads.nc chirps-v2.0.$1.dekads.nc3
gdal_translate -of netcdf -co "WRITE_BOTTOMUP=NO" chirps-v2.0.$1.dekads.nc3 chirps-v2.0.$1.dekads.flipped.nc3
nccopy -d 5 -7 -c time/4,lat/250,lon/900 chirps-v2.0.$1.dekads.flipped.nc3 chirps-v2.0.$1.dekads.nc
rm chirps-v2.0.$1.dekads.nc3
