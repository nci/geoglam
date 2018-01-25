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
#import medoid as md
import json
import os
import sys
import datetime 
import time
import argparse
# For deployment, We no longer need pyximport but pre-built fast_medoid library
#import pyximport; pyximport.install()
from fast_medoid import medoid

def get_output_filename(output_dir, source_file, ver):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    parts = source_file.split('.')
    year = parts[-3]
    hv = parts[-4]
    filename = 'FC_Monthly_Medoid.v%s.MCD43A4.%s.%s.006.nc' % (ver, hv, year)

    full_filename = os.path.join(output_dir, filename)
    return full_filename

def pack_data(src_filename, output_filename, arr, timestamps):
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

        var = dest.createVariable("phot_veg", "u1", ("time", "y", "x"), fill_value=255, zlib=True, chunksizes=(1, 240, 240))
        var.long_name = "Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = arr[:, 0, :, :]

        var = dest.createVariable("nphot_veg", "u1", ("time", "y", "x"), fill_value=255, zlib=True, chunksizes=(1, 240, 240))
        var.long_name = "Non Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = arr[:, 1, :, :]

        var = dest.createVariable("bare_soil", "u1", ("time", "y", "x"), fill_value=255, zlib=True, chunksizes=(1, 240, 240))
        var.long_name = "Bare Soil"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = arr[:, 2, :, :]

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

def compute_monthly_medoid(filename, output_dir, ver):
    with netCDF4.Dataset(filename) as ds:
        month_idx_tmp = [[] for _ in xrange(12)]    
        for its, ts in enumerate(ds['time']):
            date = netCDF4.num2date(ts, ds['time'].units)
            month_idx_tmp[date.month-1].append(its)

        #The following checks that if the current month is complete
        #We only compute medoids if the current month is complete
        for i in xrange(len(month_idx_tmp)):
          m = month_idx_tmp[i]
          if len(m) == 0:
            continue

          ts = netCDF4.num2date(ds['time'][m[-1]], ds['time'].units) 
          next_date = ts + datetime.timedelta(days=8)          
          if ts.month == next_date.month:
            month_idx_tmp[i] = []

        month_idx = [m for m in month_idx_tmp if len(m) > 0]        
        
        pv = np.asarray(ds['phot_veg'][:], dtype=np.float32)
        npv = np.asarray(ds['nphot_veg'][:], dtype=np.float32)
        soil = np.asarray(ds['bare_soil'][:], dtype=np.float32)
        
        data = np.empty((pv.shape[0], 3, pv.shape[1], pv.shape[2]), dtype=pv.dtype)
        data[:, 0, :, :] = pv
        data[:, 1, :, :] = npv
        data[:, 2, :, :] = soil

        monthly_data_list = []
        timestamp_list = []
        
        for m, m_idx in enumerate(month_idx):
            t0 = time.time()
            d = data[m_idx, ...]
            month_medoid = medoid(d)
            monthly_data_list.append(month_medoid)
            ts = netCDF4.num2date(ds['time'][m_idx[0]], ds['time'].units) 
            timestamp_list.append(ts)
            month = m + 1
            
            print '  month:%d, org_mean:%.4f, org_min:%.4f, org_max:%.4f' % (month, d.mean(), d.min(), d.max())
            print '            med_mean:%.4f, med_min:%.4f, med_max:%.4f' % (month_medoid.mean(), month_medoid.min(), month_medoid.max())
            print '            time:%.2f secs' % (time.time() - t0)

        if len(monthly_data_list) == 0:
            print 'no data for %s' % filename
            return
            
        if len(monthly_data_list) == 1:
            monthly_data = np.expand_dims(monthly_data_list[0], axis=0)
        else:
            monthly_data = np.stack(monthly_data_list)
        output_filename = get_output_filename(output_dir, filename, ver)    
        pack_data(filename, output_filename, monthly_data, timestamp_list)
    
        print '  output: %s' % output_filename
 

if __name__ == "__main__":
    #src_root_dir = '/g/data2/u39/public/prep/modis-fc/global_fc_006'
    parser = argparse.ArgumentParser(description="""Modis Vegetation Analysis argument parser""")
    
    parser.add_argument(dest="src_root_dir", type=str, help="source root path")
    parser.add_argument(dest="h", type=str, help="h of source file")
    parser.add_argument(dest="v", type=str, help="v of source file")
    parser.add_argument(dest="year", type=int, help="Year of source file")
    parser.add_argument(dest="output_dir", type=str, help="Full path to destination.")
    parser.add_argument("--version", default='310', type=str, help="Product version")
    args = parser.parse_args()

    src_root_dir = args.src_root_dir
    src_filename = os.path.join(src_root_dir, 'FC.v%s.MCD43A4.h%sv%s.%s.006.nc' % (args.version, args.h, args.v, args.year))
    
    print '%s' % src_filename
    compute_monthly_medoid(src_filename, args.output_dir, args.version)
    

