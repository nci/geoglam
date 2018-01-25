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
import json
import os
import sys
import datetime 
import time
import argparse

def get_output_filename(output_dir, source_file, year, product_name, ver):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    parts = source_file.split('.')
    hv = parts[-4]
    filename = 'FC_%s.v%s.MCD43A4.%s.%s.006.nc' % (product_name, ver, hv, year)

    full_filename = os.path.join(output_dir, filename)
    return full_filename

def pack_data(src_filename, output_filename, pv_pct, npv_pct, soil_pct, data_type, timestamps):
    with netCDF4.Dataset(output_filename, 'w', format='NETCDF4') as dest:
        with open('nc_metadata.json') as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(dest, key, attrs[key])
        setattr(dest, "date_created", datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
    
        ds = gdal.Open('NETCDF:"%s":phot_veg' % src_filename)
        proj_wkt = ds.GetProjection()
        geot = ds.GetGeoTransform()

        t_dim = dest.createDimension("time", len(timestamps))
        x_dim = dest.createDimension("x", ds.RasterXSize)
        y_dim = dest.createDimension("y", ds.RasterYSize)

        var = dest.createVariable("time", "f8", ("time",))
        var.units = "seconds since 1970-01-01 00:00:00.0"
        var.calendar = "standard"
        var.long_name = "Time, unix time-stamp"
        var.standard_name = "time"
        var[:] = netCDF4.date2num(timestamps, units="seconds since 1970-01-01 00:00:00.0", calendar="standard")

        var = dest.createVariable("x", "f8", ("x",))
        var.units = "m"
        var.long_name = "x coordinate of projection"
        var.standard_name = "projection_x_coordinate"
        var[:] = np.linspace(geot[0], geot[0]+(geot[1]*ds.RasterXSize), ds.RasterXSize)

        var = dest.createVariable("y", "f8", ("y",))
        var.units = "m"
        var.long_name = "y coordinate of projection"
        var.standard_name = "projection_y_coordinate"
        var[:] = np.linspace(geot[3], geot[3]+(geot[5]*ds.RasterYSize), ds.RasterYSize)

        var = dest.createVariable("phot_veg", data_type, ("time", "y", "x"), fill_value=255, zlib=True, chunksizes=(1, 240, 240))
        var.long_name = "Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = pv_pct

        var = dest.createVariable("nphot_veg", data_type, ("time", "y", "x"), fill_value=255, zlib=True, chunksizes=(1, 240, 240))
        var.long_name = "Non Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = npv_pct

        var = dest.createVariable("bare_soil", data_type, ("time", "y", "x"), fill_value=255, zlib=True, chunksizes=(1, 240, 240))
        var.long_name = "Bare Soil"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = soil_pct

        var = dest.createVariable("sinusoidal", 'S1', ())

        var.grid_mapping_name = "sinusoidal"
        var.false_easting = 0.0
        var.false_northing = 0.0
        var.longitude_of_central_meridian = 0.0
        var.longitude_of_prime_meridian = 0.0
        var.semi_major_axis = 6371007.181
        var.inverse_flattening = 0.0
        var.spatial_ref = proj_wkt
        var.GeoTransform = "{} {} {} {} {} {} ".format(*[geot[i] for i in range(6)])

def combine_by_month(src_root_filename, input_dir, src_file, h, v, year_start, year_end, output_dir, product_name, data_type, ver):
    for i_year in xrange(year_start, year_end+1):
        print 'processing year: %d' % i_year
        pv_data = npv_data = soil_data = None
        timestamp_list = []
        t0 = time.time()
        for i_month in xrange(0, 12):
            month = i_month + 1
            year = i_year
            src_filename = os.path.join(input_dir, src_file % (product_name, h, v, month, year))
    
            if not os.path.isfile(src_filename):
                print 'monthly percentile file not found: "%s"' % src_filename
                continue

            with netCDF4.Dataset(src_filename) as ds:
                if data_type == 'f4':
                    dtype = np.float32
                elif data_type == 'u1':
                    dtype = np.uint8
                else:
                    raise Exception('invalid data type')
                pv = np.asarray(ds['phot_veg'][:], dtype=dtype)
                npv = np.asarray(ds['nphot_veg'][:], dtype=dtype)
                soil = np.asarray(ds['bare_soil'][:], dtype=dtype)

                if pv_data is None:
                    pv_data = np.zeros((12, pv.shape[1], pv.shape[2]), dtype=pv.dtype)
                    npv_data = np.zeros_like(pv_data)
                    soil_data = np.zeros_like(pv_data)

                pv_data[i_month, ...] = np.squeeze(pv)
                npv_data[i_month, ...] = np.squeeze(npv)
                soil_data[i_month, ...] = np.squeeze(soil)

                ts = netCDF4.num2date(ds['time'][0], ds['time'].units) 
                timestamp_list.append(ts)

        if pv_data is None:
            print 'no data, month:%d, src_root:%s' % (i_month, src_root_filename)
            continue

        print time.time() - t0
        
        output_filename = get_output_filename(output_dir, src_root_filename, year, product_name, ver)
        print output_filename
        pack_data(src_root_filename, output_filename, pv_data, npv_data, soil_data, data_type, timestamp_list)

if __name__ == "__main__":
    #src_root_dir = '/g/data2/u39/public/prep/modis-fc/FC.v302.MCD43A4'
    parser = argparse.ArgumentParser(description="""Modis Vegetation Analysis argument parser""")
    parser.add_argument(dest="src_root_dir", type=str, help="Source root dir")
    parser.add_argument(dest="input_dir", type=str, help="Full path of input file")    
    parser.add_argument(dest="h", type=str, help="h of source file")
    parser.add_argument(dest="v", type=str, help="v of source file")
    parser.add_argument(dest="year_start", type=int, help="Starting year of source file")
    parser.add_argument(dest="year_end", type=int, help="Ending year of source file (inclusive)")
    parser.add_argument(dest="output_dir", type=str, help="Full path to destination.")
    parser.add_argument(dest="product_name", type=str, help="Product name")
    parser.add_argument(dest="data_type", type=str, help="Product data type")
    parser.add_argument("--version", default='310', type=str, help="Product version")
    args = parser.parse_args()

    src_root_dir = args.src_root_dir
    src_file = 'FC_%s.v' + args.version + '.MCD43A4.h%sv%s.%s.%s.006.nc'
    
    src_root_filename = os.path.join(src_root_dir, 'FC.v%s.MCD43A4.h%sv%s.%s.006.nc' % (args.version, args.h, args.v, args.year_start))
    combine_by_month(src_root_filename, args.input_dir, src_file, args.h, args.v, args.year_start, args.year_end, args.output_dir, args.product_name, args.data_type, args.version)
    

