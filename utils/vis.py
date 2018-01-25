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
import numpy as np
from osgeo import gdal
import netCDF4
import matplotlib.pyplot as plt
import argparse

if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description="""Modis Vegetation Analysis argument parser""")
    parser.add_argument(dest="source", type=str, help="Full path to source.")
    args = parser.parse_args()

    src_path = args.source
    
    f = gdal.Open(src_path)

    datasets = f.GetSubDatasets()

    band_list = []
    for ds in datasets:
        data = gdal.Open(ds[0]).ReadAsArray()
        print data[0, :, :].mean(), data[0, :, :].max(), data[0, :, :].min()
        print data.shape, data.dtype
        band_list.append(data[0, :, :])

    data = np.dstack(band_list)

    plt.imshow(data)
    plt.show()

