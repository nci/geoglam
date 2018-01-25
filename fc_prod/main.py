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
import scipy.optimize as opt
import scipy.ndimage
from PIL import Image
import netCDF4
import json
import os
import sys
import datetime
import argparse
import time

def get_nc_name(hdf_tile_path, dest, ver, suffix=""):
    file_name = hdf_tile_path.split("/")[-1][:-4]
    file_name_parts = file_name.split(".")
    folder =  ".".join([file_name_parts[2], file_name_parts[1][1:-3]])
    dest_fold = os.path.join(dest, folder)
    if not os.path.exists(dest_fold):
        os.makedirs(dest_fold)
    short_file_name = file_name.split(".")[:-1]
    return os.path.join(dest_fold, "FC" + suffix + (".v%s." % ver) + ".".join(short_file_name) + ".nc")


#def generate_thumbnail(arr, file_name):
#    norm_arr = (arr.astype(np.float32).clip(0,100)*2.55).astype(np.uint8)
#    img = Image.fromarray(norm_arr, 'RGB')
#    g, b, r = img.split()
#    img = Image.merge("RGB", (r, g, b))
#    img.save(file_name[:-2] + "png")


def pack_data_lowres(hdf_file, arr, dest, ver):
    file_name = get_nc_name(hdf_file, dest, ver, "_LR")
    
    with netCDF4.Dataset(file_name, 'w', format='NETCDF4') as dest:
        with open('nc_metadata.json') as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(dest, key, attrs[key])
        setattr(dest, "date_created", datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
        
        ds = gdal.Open('HDF4_EOS:EOS_GRID:"{}":MOD_Grid_BRDF:Nadir_Reflectance_Band1'.format(hdf_file))
        proj_wkt = ds.GetProjection()
        geot = ds.GetGeoTransform()
        geot_p = ds.GetGeoTransform()
        geot = (geot_p[0], geot_p[1]*10, geot_p[2], geot_p[3], geot_p[4], geot_p[5]*10)
        
        x_dim = dest.createDimension("x", ds.RasterXSize/10)
        y_dim = dest.createDimension("y", ds.RasterYSize/10)

        var = dest.createVariable("x", "f8", ("x",))
        var.units = "m"
        var.long_name = "x coordinate of projection"
        var.standard_name = "projection_x_coordinate"
        var[:] = np.linspace(geot[0], geot[0]+(geot[1]*ds.RasterXSize), ds.RasterXSize/10)
        
        var = dest.createVariable("y", "f8", ("y",))
        var.units = "m"
        var.long_name = "y coordinate of projection"
        var.standard_name = "projection_y_coordinate"
        var[:] = np.linspace(geot[3], geot[3]+(geot[5]*ds.RasterYSize), ds.RasterYSize/10)
        
        var = dest.createVariable("b", "f8", ("b",))
        var.units = "%"
        var.long_name = "photosynthetic vegetation, non-photosynthetic vegetation, bare soil"
        var.standard_name = "bands"
        #var[:] = np.array(["photosynthetic vegetation", "non-photosynthetic vegetation", "bare soil"])
        var[:] = np.array([1, 2, 3])
        
        resampled = scipy.ndimage.zoom(arr, [.1, .1, 1])
 
        var = dest.createVariable("phot_veg", "u1", ("y", "x"), fill_value=255)
        var.long_name = "Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = resampled[:, :, 0]
        
        var = dest.createVariable("nphot_veg", "u1", ("y", "x"), fill_value=255)
        var.long_name = "Non Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = resampled[:, :, 1]
       
        var = dest.createVariable("bare_soil", "u1", ("y", "x"), fill_value=255)
        var.long_name = "Bare Soil"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = resampled[:, :, 2]
        
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

def pack_data(hdf_file, arr, dest, ver):
    
    file_name = get_nc_name(hdf_file, dest, ver)
    
    with netCDF4.Dataset(file_name, 'w', format='NETCDF4') as dest:
        with open('nc_metadata.json') as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(dest, key, attrs[key])
        setattr(dest, "date_created", datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
        
        ds = gdal.Open('HDF4_EOS:EOS_GRID:"{}":MOD_Grid_BRDF:Nadir_Reflectance_Band1'.format(hdf_file))
        proj_wkt = ds.GetProjection()
        geot = ds.GetGeoTransform()
        
        x_dim = dest.createDimension("x", ds.RasterXSize)
        y_dim = dest.createDimension("y", ds.RasterYSize)

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
        
        var = dest.createVariable("phot_veg", "u1", ("y", "x"), fill_value=255)
        var.long_name = "Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = arr[:, :, 0]
        
        var = dest.createVariable("nphot_veg", "u1", ("y", "x"), fill_value=255)
        var.long_name = "Non Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = arr[:, :, 1]
       
        var = dest.createVariable("bare_soil", "u1", ("y", "x"), fill_value=255)
        var.long_name = "Bare Soil"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = arr[:, :, 2]
        
        #generate_thumbnail(arr, file_name)
        
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

def jp_superfunc(A, arr, mask):
    res = np.zeros((arr.shape[0], A.shape[1]), dtype=np.uint8)
    
    for i in range(arr.shape[0]):
        #if mask[i]:
        if mask[i] and (arr[i, :7] < 0).sum() == 0 and (arr[i, :7] > 1).sum() == 0:
            res[i, :] = (opt.nnls(A, arr[i, :])[0].clip(0, 2.54)*100).astype(np.uint8)
        else:
            res[i, :] = np.ones((A.shape[1]), dtype=np.uint8)*(255) 

    return res

def input_mask(pq_mask_path):
    #pq_ds = gdal.Open('HDF4_EOS:EOS_GRID:"{}":MOD_Grid_BRDF:BRDF_Albedo_Ancillary'.format(pq_mask_path))
    pq_ds = gdal.Open('HDF4_EOS:EOS_GRID:"{}":MOD_Grid_BRDF:BRDF_Albedo_LandWaterType'.format(pq_mask_path))
    pq_raw = pq_ds.GetRasterBand(1).ReadAsArray()

    mask = 0b0000000000001111
    #pq = pq_raw>>4 & mask
    # in collection 6 I don't need to shift 4 places  
    pq = pq_raw & mask

    snow_ds = gdal.Open('HDF4_EOS:EOS_GRID:"{}":MOD_Grid_BRDF:Snow_BRDF_Albedo'.format(pq_mask_path))
    snow_pq = np.equal(snow_ds.ReadAsArray(), 0)

    # True if pq is 1, 2 or 4    
    pq = np.logical_and(snow_pq, np.logical_or(np.logical_or(np.equal(pq, np.ones(1)), np.equal(pq, np.ones(1)*2)), np.equal(pq, np.ones(1)*4)))
    #pq = np.logical_or(np.logical_or(np.equal(pq, np.ones(1)), np.equal(pq, np.ones(1)*2)), np.equal(pq, np.ones(1)*4))
    
    return pq.reshape(2400*2400)

def input_stack(tile_path):
    bands = [gdal.Open('HDF4_EOS:EOS_GRID:"{}":MOD_Grid_BRDF:Nadir_Reflectance_Band{}'.format(tile_path, b)) for b in range(1, 8)]

    bands_stack = np.empty((2400*2400, 85), dtype=np.float32)
    
    #add band
    for b in range(0, len(bands)):
        bands_stack[:, b] = np.nan_to_num( bands[b].GetRasterBand(1).ReadAsArray().reshape((2400*2400, ))*.0001 )

    ii = 7
    #add band*log(band)
    bands_stack[:, ii:ii+7] = np.nan_to_num( np.log(bands_stack[:, :7]) )
    
    #add band*log(band)
    bands_stack[:, ii+7:ii+14] = np.nan_to_num( bands_stack[:, :7] * bands_stack[:, ii:ii+7] )

    ii += 14
    #add band*next_band
    for b in range(7):
        for b2 in range(b+1, 7):
            bands_stack[:, ii] = np.nan_to_num( bands_stack[:, b] * bands_stack[:, b2] )
            ii += 1
    
    #add log(band)*log(next_band)
    for b in range(7):
        for b2 in range(b+1, 7):
            bands_stack[:, ii] = np.nan_to_num( bands_stack[:, b+7] * bands_stack[:, b2+7] )
            ii += 1
    
    #add (next_band-band)/(next_band+band)
    for b in range(7):
        for b2 in range(b+1, 7):
            bands_stack[:, ii] = np.nan_to_num( (bands_stack[:, b2]-bands_stack[:, b]) / (bands_stack[:, b2]+bands_stack[:, b]) )
            ii += 1
    
    #--------------------------------
    # this bit new, by JPG
    bands_stack[:, -1] = np.ones((bands_stack.shape[0],), dtype=bands_stack.dtype)
    #--------------------------------
    
    return bands_stack


def members():

    #A = np.load("endmembers.npy")
    A = np.loadtxt("endmembers_v6_20170831.txt")
    #--------------------------------
    # this bit new, by JPG
    SumToOneWeight = 0.02
    ones = np.ones(A.shape[1]) * SumToOneWeight
    ones = ones.reshape(1, A.shape[1])
    A = np.concatenate((A, ones), axis=0).astype(np.float32)

    return A

 
if __name__ == "__main__":
    
    #root_path = "/g/data2/u39/public/data/modis/lpdaac-tiles-c6/"

    parser = argparse.ArgumentParser(description="""Modis Vegetation Analysis argument parser""")
    parser.add_argument(dest="root_path", type=str, help="Root path")    
    parser.add_argument(dest="year", type=int, help="Year of a modis tile (HDF-EOS).")
    parser.add_argument(dest="month", type=int, help="Year of a modis tile (HDF-EOS).")
    #parser.add_argument(dest="julday", type=int, help="Julian Day of a modis tile (HDF-EOS).")
    parser.add_argument(dest="h", type=str, help="Horizontal Corrrdinate of Modis Tile (HDF-EOS).")
    parser.add_argument(dest="v", type=str, help="Vertical Corrrdinate of Modis Tile (HDF-EOS).")
    parser.add_argument(dest="destination", type=str, help="Full path to destination.")
    parser.add_argument("--version", default='310', type=str, help="Product version")
    parser.add_argument("--update_file", default='', type=str, help="Existing prodcut file to update")
    args = parser.parse_args()

    tile_year = args.year
    tile_month = args.month
    root_path = args.root_path
    #tile_julday = args.julday
    tile_h = args.h
    tile_v = args.v
    dest_path = args.destination

    ts_max = None
    if os.path.isfile(args.update_file):
        print 'Updating %s' % args.update_file
        with netCDF4.Dataset(args.update_file, 'r', format='NETCDF4') as src:
            ts = netCDF4.num2date(src['time'][:], src['time'].units, src['time'].calendar)
            ts_max = np.datetime64(ts.max())
            print 'latest existing timestamp to update: %s' % str(ts_max)
            
            tile_ts = np.datetime64('%d-%02d' % (tile_year, tile_month), 'D')
            print 'tile timestamp: %s' % str(tile_ts)
            
            if tile_ts <= ts_max:
                print 'tile timestamp is older than existing timestamp of the product. No updates needed.'
                sys.exit(0)
    else:
      print 'Update file not found: %s, creating FC products from scratch' % args.update_file

    data_path = os.path.join(root_path, "MCD43A4.006")
    tile_path = [os.path.join(data_path, p) for p in os.listdir(data_path) if p.startswith('%d.%.2d' % (tile_year, tile_month))]
        
    days = []
    for t, tf in enumerate(tile_path):
        timestamp = tf.split('/')[-1]
        ts_parts = timestamp.split('.')
        day = int(ts_parts[-1])
        
        if ts_max is not None:
            year = int(ts_parts[0])
            month = int(ts_parts[1])
            
            ts = np.datetime64('%d-%02d-%02d' % (year, month, day))
            if ts <= ts_max:
                continue
                
        days.append(day)
    
    #we cannot assume the first day of the month starts from 1. For example,
    #Feb in 2016 starts from 2 in this dataset
    #Thus we need to manually extract the first day of the month and then extract every 8th day from there     
    sorted_days = np.sort(days)
    desired_days = []
    last_day = -1
    for i, sd in enumerate(sorted_days):
        if i == 0:
            desired_days.append(sd)
            last_day = sd
        else:
            if sd - last_day >= 8:
                desired_days.append(sd)
                last_day = sd

    desired_timestamps = []    
    for d in desired_days:
        ts = '%d.%.2d.%.2d' % (tile_year, tile_month, d)  
        desired_timestamps.append(ts)
    
    # We will need to overlay the A2 mask on the A4 data
    # So we first loop through the timestamps then filter the data files
    # to ensure of the order of data files
    _paths = [os.path.join(data_path, p) for p in os.listdir(data_path)]
    desired_tile_path = []
    for dt in desired_timestamps:
        for p in _paths:
            if dt in p:
                desired_tile_path.append(p)
                break
    
    _tile_files = []
    for tp in desired_tile_path:
        tile_found = False
        for t in os.listdir(tp):
            if '.h%sv%s.' % (tile_h, tile_v) in t and t.endswith('.hdf'):
                tf = os.path.join(tp, t)
                _tile_files.append(tf)
                tile_found = True
                break
        if not tile_found:
          print 'A4 tile h%sv%s data file not found at %s' % (tile_h, tile_v, tp)
          _tile_files.append(None)

    
    mask_data_path = os.path.join(root_path, "MCD43A2.006")
    _paths = [os.path.join(mask_data_path, p) for p in os.listdir(mask_data_path)]
    mask_tile_path = []
    for dt in desired_timestamps:
        for p in _paths:
            if dt in p:
                mask_tile_path.append(p)
                break
                
    _mask_tile_files = []
    for tp in mask_tile_path:
        tile_found = False
        for t in os.listdir(tp):
            if '.h%sv%s.' % (tile_h, tile_v) in t and t.endswith('.hdf'):
                tf = os.path.join(tp, t)
                _mask_tile_files.append(tf)
                tile_found = True
                break
        if not tile_found:
          print 'A2 tile h%sv%s data file not found at %s' % (tile_h, tile_v, tp)
          _mask_tile_files.append(None)
      
    assert len(_tile_files) == len(_mask_tile_files), \
      'tile_files, mask_tile_files length not equal. %d, %d, %s, %s' % (args.year, 
      args.month, str(_tile_files), str(_mask_tile_files))

    tile_files = []
    mask_tile_files = []

    for im, modis_tile in enumerate(_tile_files):
      if modis_tile is not None and _mask_tile_files[im] is not None:
        tile_files.append(modis_tile)
        mask_tile_files.append(_mask_tile_files[im])

    A = members()
    for im, modis_tile in enumerate(tile_files):
        print 'processing tile: %s (%d of %d):' % (modis_tile, im+1, len(tile_files))
        print 'processing mask tile:%s (%d of %d):' % (mask_tile_files[im], im+1, len(tile_files))
        t0 = time.time()
        bands_stack = input_stack(modis_tile)
        
        print '    time input_stack(): ', time.time() - t0
           
        mask = input_mask(mask_tile_files[im])
        
        t0 = time.time()
        out = jp_superfunc(A, bands_stack, mask).reshape(2400, 2400, 3)
        print '    time jp_superfunc(): ', time.time() - t0
        
        print 'saving output to %s' % dest_path
        pack_data(modis_tile, out, dest_path, args.version)
 

    
