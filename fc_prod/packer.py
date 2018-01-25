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
import functools
import glob
import json

import netCDF4
import numpy as np
import os.path

def pack(infoldername, outfoldername, ver, low_res=False, update_file=''):
    
    id = os.path.basename(os.path.normpath(infoldername))
    
    # Config - Where are the source files?  What resolution to use?
    config = {
        True: {  # low_res=True, use lowres tiles
            "source_path": os.path.join(infoldername + "/FC_LR.*.nc"),
            "dest_filename": os.path.join(outfoldername, "FC_LR.v{ver}.MCD43A4.{id}.006.nc"),
            "scale": 10.0},
        False: {  # low_res=False, use full resolution tiles
            "source_path": os.path.join(infoldername + "/FC.*.nc"),
            "dest_filename": os.path.join(outfoldername, "FC.v{ver}.MCD43A4.{id}.006.nc"),
            "scale": 1.0}
        }

    timestamps = []
    source_files = glob.glob(config[low_res]["source_path"].format(id))
    if len(source_files) == 0:
        print 'No source files at %s' % config[low_res]["source_path"].format(id)
        return 
        
    if not os.path.isfile(update_file):
        phot_stack = None
        nphot_stack = None
        bare_stack = None
        
        print 'Update file not found: %s, creating FC products from scratch' % update_file
    else:
        with netCDF4.Dataset(update_file, 'r', format='NETCDF4') as src:
            phot_stack = src["phot_veg"][:]
            nphot_stack = src["nphot_veg"][:]
            bare_stack = src["bare_soil"][:]
            ts = netCDF4.num2date(src['time'][:], src['time'].units, src['time'].calendar)
            for i in xrange(ts.shape[0]):
                timestamps.append(ts[i])

        print 'Updating %s' % update_file, phot_stack.shape    
        
    proj_wkt = None
    geot = None

    for file in sorted(source_files):
        tile_ts = file.split("/")[-1].split(".")[3][1:]
        date = datetime.datetime(int(tile_ts[:4]), 1, 1) + datetime.timedelta(int(tile_ts[4:])-1)
        timestamps.append(date)
        
        with netCDF4.Dataset(file, 'r', format='NETCDF4') as src:
            if geot is None:
                var = src["sinusoidal"]
                proj_wkt = var.spatial_ref
                geot = [float(val) for val in var.GeoTransform.split(" ") if val]

            if phot_stack is None:
                phot_stack = np.expand_dims(src["phot_veg"][:], axis=0)
                nphot_stack = np.expand_dims(src["nphot_veg"][:], axis=0)
                bare_stack = np.expand_dims(src["bare_soil"][:], axis=0)
            else:
                phot_stack = np.vstack((phot_stack, np.expand_dims(src["phot_veg"][:], axis=0)))
                nphot_stack = np.vstack((nphot_stack, np.expand_dims(src["nphot_veg"][:], axis=0)))
                bare_stack = np.vstack((bare_stack, np.expand_dims(src["bare_soil"][:], axis=0)))

    assert phot_stack is not None
    print 'Saving file %s' % config[low_res]["dest_filename"].format(ver=ver,id=id), phot_stack.shape
    with netCDF4.Dataset(config[low_res]["dest_filename"].format(ver=ver, id=id), 'w', format='NETCDF4') as dest:
        with open('nc_metadata.json') as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(dest, key, attrs[key])
 
        setattr(dest, "date_created", datetime.datetime.utcnow().isoformat())
        
        t_dim = dest.createDimension("time", len(timestamps))
        x_dim = dest.createDimension("x", phot_stack.shape[2])
        y_dim = dest.createDimension("y", phot_stack.shape[1])
        
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
        var[:] = np.linspace(
            geot[0],
            geot[0] + (config[low_res]["scale"] * geot[1] * phot_stack.shape[2]),
            phot_stack.shape[2])

        var = dest.createVariable("y", "f8", ("y",))
        var.units = "m"
        var.long_name = "y coordinate of projection"
        var.standard_name = "projection_y_coordinate"
        var[:] = np.linspace(geot[3], geot[3]+(config[low_res]["scale"]*geot[5]*phot_stack.shape[1]), phot_stack.shape[1])

        # Same variable properties for all data.
        # NB: use chunksizes = (5, 240, 240) for low-res data.
        create_data_var = functools.partial(
            dest.createVariable, datatype='u1', dimensions=("time", "y", "x"),
            fill_value=255, zlib=True, chunksizes=(1, 240, 240))

        var = create_data_var("phot_veg")
        var.long_name = "Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = phot_stack

        var = create_data_var("nphot_veg")
        var.long_name = "Non Photosynthetic Vegetation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = nphot_stack

        var = create_data_var("bare_soil")
        var.long_name = "Bare Soil"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = bare_stack

        var = dest.createVariable("sinusoidal", 'S1', ())
        var.grid_mapping_name = "sinusoidal"
        var.false_easting = 0.0
        var.false_northing = 0.0
        var.longitude_of_central_meridian = 0.0
        var.longitude_of_prime_meridian = 0.0
        var.semi_major_axis = 6371007.181
        var.inverse_flattening = 0.0
        var.spatial_ref = proj_wkt
        var.GeoTransform = "{} {} {} {} {} {}".format(geot[0], config[low_res]["scale"]*geot[1], geot[2], geot[3], geot[4], config[low_res]["scale"]*geot[5])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Modis Vegetation Analysis NetCDF aggregator.")
    parser.add_argument(dest="infoldername", type=str, help="Input folder to pack.")
    parser.add_argument(dest="outfoldername", type=str, help="Output folder to write.")
    parser.add_argument("--version", default='310', type=str, help="Product version")
    parser.add_argument("--update_file", default='', type=str, help="Existing prodcut file to update")
    args = parser.parse_args()

    infoldername = args.infoldername
    outfoldername = args.outfoldername

    pack(infoldername, outfoldername, args.version, False, update_file=args.update_file)
    #pack(infoldername, outfoldername, args.version, True, update_file=args.update_file)
