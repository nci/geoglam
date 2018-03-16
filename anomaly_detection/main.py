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

def get_output_filename(output_dir, source_file, month, product_name, ver):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    parts = source_file.split('.')
    year = parts[-3]
    hv = parts[-4]
    filename = 'FC_%s.v%s.MCD43A4.%s.%d.%s.006.nc' % (product_name, ver, hv, month, year)

    full_filename = os.path.join(output_dir, filename)
    return full_filename

def pack_data(src_filename, output_filename, pv_data, npv_data, soil_data, data_type, timestamps):
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

        var = dest.createVariable("phot_veg", data_type, ("time", "y", "x"), fill_value=255, zlib=True)
        var.long_name = "Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = pv_data

        var = dest.createVariable("nphot_veg", data_type, ("time", "y", "x"), fill_value=255, zlib=True)
        var.long_name = "Non Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = npv_data

        var = dest.createVariable("bare_soil", data_type, ("time", "y", "x"), fill_value=255, zlib=True)
        var.long_name = "Bare Soil"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = soil_data

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

def compute_percentiles(data):
    sort_idx = np.zeros_like(data)
    idx = np.arange(sort_idx.shape[0])
    
    arg_idx = np.argsort(data, axis=0)
    for i in xrange(sort_idx.shape[1]):
        for j in xrange(sort_idx.shape[2]):
           sort_idx[arg_idx[:, i, j], i, j] = idx 

    masks = data < 255
    count = masks.sum(axis=0)
    #The original formula is percentile = ((p-1)/(n-1)) * 100 
    #But numpy index starts from zero. So we use p/(n-1) * 100 instead
    percentiles = sort_idx / (count - 1.) * 100
    percentiles[~masks] = 255
    nan_masks = np.isnan(percentiles) | np.isinf(percentiles)
    percentiles[nan_masks] = 255
    percentiles[percentiles < 0] = 255
    
    #print percentiles.max(), percentiles.min(), percentiles.shape
    #print sort_idx.max(), sort_idx.min(), count.max(), count.min(), percentiles.max(), percentiles.min()
    percentiles = np.round(percentiles).astype(np.uint8)
    return percentiles

def compute_mean_diff(data):
    masks = data < 255
    _data = data.copy()
    _data[~masks] = 0
    data_sum = _data.sum(axis=0)
    counts = masks.sum(axis=0)
    mean_data = data_sum / counts
    diff = data - mean_data
    nan_masks = np.isnan(diff) | np.isinf(diff)
    diff[nan_masks] = 255
    diff[~masks] = 255
    
    return diff
                
def compute_by_month(src_root_filename, raw_src_file_list, month, mean_diffs_output_dir, percentiles_output_dir, ver):
    pv_data = npv_data = soil_data = None
    t0 = time.time()
    
    src_file_list = []
    for src_f in raw_src_file_list:
      if os.path.isfile(src_f):
        src_file_list.append(src_f)
      else:
        print 'source file not found: %s' % src_f
    
    timestamp_list = []
    for _, src_f in enumerate(src_file_list):
      parts = src_f.split('.')
      year = int(parts[-3])
      timestamp_list.append(datetime.datetime(year, month, 1, 0, 0))

    for i_src, src_f in enumerate(src_file_list):
        with netCDF4.Dataset(src_f) as ds:
            month2idx = -1
            for its, ts in enumerate(ds['time']):
                date = netCDF4.num2date(ts, ds['time'].units)
                if date.month == month:
                    month2idx = its
                    break

            if month2idx == -1:
                print '%s has no data for month: %d' % (src_f, month)
                continue            
            
            pv = np.asarray(ds['phot_veg'][month2idx, ...], dtype=np.float32)
            npv = np.asarray(ds['nphot_veg'][month2idx, ...], dtype=np.float32)
            soil = np.asarray(ds['bare_soil'][month2idx, ...], dtype=np.float32)
            ts = netCDF4.num2date(ds['time'][month2idx], ds['time'].units) 

            if pv_data is None:
                fill_val = 255
                pv_data = fill_val * np.ones((len(src_file_list), pv.shape[0], pv.shape[1]), dtype=pv.dtype)
                npv_data = fill_val * np.ones_like(pv_data)
                soil_data = fill_val * np.ones_like(pv_data)

            pv_data[i_src, ...] = pv
            npv_data[i_src, ...] = npv
            soil_data[i_src, ...] = soil
            timestamp_list[i_src] = ts

    if pv_data is None:
        raise Exception('no data, month:%d, src_root:%s' % (month, src_root_filename))

    pv_mean_diff = compute_mean_diff(pv_data)
    npv_mean_diff = compute_mean_diff(npv_data)
    soil_mean_diff = compute_mean_diff(soil_data)
    
    #we save one file per year for the given month. This simplifies combine_outputs.py
    for i_src, src_f in enumerate(src_file_list):
        output_filename = get_output_filename(mean_diffs_output_dir, src_f, month, 'Mean_Diff', ver)
        print output_filename
        #print pv_ranks[i_src, ...].shape
        pack_data(src_root_filename, 
            output_filename, 
            np.expand_dims(pv_mean_diff[i_src, ...], 0), 
            np.expand_dims(npv_mean_diff[i_src, ...], 0), 
            np.expand_dims(soil_mean_diff[i_src, ...], 0), 'f4', [timestamp_list[i_src], ])

    pv_pct = compute_percentiles(pv_data)
    npv_pct = compute_percentiles(npv_data)
    soil_pct = compute_percentiles(soil_data)
   
    for i_src, src_f in enumerate(src_file_list):
        output_filename = get_output_filename(percentiles_output_dir, src_f, month, 'Percentile', ver)
        print output_filename
        #print pv_ranks[i_src, ...].shape
        pack_data(src_root_filename, 
            output_filename, 
            np.expand_dims(pv_pct[i_src, ...], 0), 
            np.expand_dims(npv_pct[i_src, ...], 0), 
            np.expand_dims(soil_pct[i_src, ...], 0), 'u1', [timestamp_list[i_src], ])
        
     
    print 'time elapsed: ', time.time() - t0


if __name__ == "__main__":
    #src_root_dir = 'src_root_dir = '/g/data2/u39/public/prep/modis-fc/global_fc_006'
    parser = argparse.ArgumentParser(description="""Modis Vegetation Analysis argument parser""")
    parser.add_argument(dest="src_root_dir", type=str, help="Source root dir")
    parser.add_argument(dest="input_dir", type=str, help="Full path of input file")    
    parser.add_argument(dest="h", type=str, help="h of source file")
    parser.add_argument(dest="v", type=str, help="v of source file")
    parser.add_argument(dest="year_start", type=int, help="Starting year of source file")
    parser.add_argument(dest="year_end", type=int, help="Ending year of source file (inclusive)")
    parser.add_argument(dest="month", type=int, help="Month of source file")
    parser.add_argument(dest="mean_diffs_output_dir", type=str, help="Full path to destination of mean differences.")
    parser.add_argument(dest="percentiles_output_dir", type=str, help="Full path to destination of percentiles.")
    parser.add_argument("--version", default='310', type=str, help="Product version")
    args = parser.parse_args()

    src_root_dir = args.src_root_dir
    src_file_list = [os.path.join(args.input_dir, 'FC_Monthly_Medoid.v%s.MCD43A4.h%sv%s.%s.006.nc' % (args.version, args.h, args.v, year)) for year in xrange(args.year_start, args.year_end+1)]

    src_root_filename = os.path.join(src_root_dir, 'FC.v%s.MCD43A4.h%sv%s.%s.006.nc' % (args.version, args.h, args.v, args.year_start))
    compute_by_month(src_root_filename, src_file_list, args.month, args.mean_diffs_output_dir, args.percentiles_output_dir, args.version)
    

