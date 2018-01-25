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
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Modis Vegetation Analysis argument parser""")
    parser.add_argument(dest="input_path", type=str, help="FC product path for the files to update")    
    parser.add_argument(dest="year", type=int, help="Year of a modis tile (HDF-EOS).")
    parser.add_argument(dest="h", type=str, help="Horizontal Corrrdinate of Modis Tile (HDF-EOS).")
    parser.add_argument(dest="v", type=str, help="Vertical Corrrdinate of Modis Tile (HDF-EOS).")
    parser.add_argument("--version", default='310', type=str, help="Product version")
    args = parser.parse_args()
    
    update_filename = 'FC.v{}.MCD43A4.h{}v{}.{}.006.nc'.format(args.version, args.h, args.v, args.year)
    update_file = os.path.join(args.input_path, update_filename)
    
    print update_file
