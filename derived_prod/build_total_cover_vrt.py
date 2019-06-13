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
import xml.etree.ElementTree as ET
from osgeo import gdal

def get_total_cover_band(idx, data_file):
  return '''
  <VRTRasterBand dataType="Byte" band="{0}" subClass="VRTDerivedRasterBand">
    <Metadata>
      <MDI key="grid_mapping">sinusoidal</MDI>
      <MDI key="long_name">Total_Cover</MDI>
      <MDI key="NETCDF_DIM_time">1514764800</MDI>
      <MDI key="NETCDF_VARNAME">total_cover</MDI>
      <MDI key="units">%</MDI>
      <MDI key="_FillValue">255</MDI>
    </Metadata>
    <NoDataValue>255</NoDataValue>
    <UnitType>%</UnitType>
    <PixelFunctionType>sum</PixelFunctionType>
    <SimpleSource>
      <SourceFilename relativeToVRT="0">NETCDF:{1}:phot_veg</SourceFilename>
      <SourceBand>{0}</SourceBand>
      <SourceProperties RasterXSize="2400" RasterYSize="2400" DataType="Byte" BlockXSize="240" BlockYSize="240" />
      <SrcRect xOff="0" yOff="0" xSize="2400" ySize="2400" />
      <DstRect xOff="0" yOff="0" xSize="2400" ySize="2400" />
    </SimpleSource>
    <SimpleSource>
      <SourceFilename relativeToVRT="0">NETCDF:{1}:nphot_veg</SourceFilename>
      <SourceBand>{0}</SourceBand>
      <SourceProperties RasterXSize="2400" RasterYSize="2400" DataType="Byte" BlockXSize="240" BlockYSize="240" />
      <SrcRect xOff="0" yOff="0" xSize="2400" ySize="2400" />
      <DstRect xOff="0" yOff="0" xSize="2400" ySize="2400" />
    </SimpleSource>
  </VRTRasterBand>
'''.format(idx, data_file)

def get_vrt_template(src_ds_file):
  src_ds = gdal.Open(args.src_ds)

  vrt_tpl_file = '/vsimem/vrt_template.vrt'
  gdal.Translate(vrt_tpl_file, src_ds)
  src_ds = None

  mem_f = gdal.VSIFOpenL(vrt_tpl_file, 'rb')

  vrt_template = ''
  while not gdal.VSIFEofL(mem_f):
    data = gdal.VSIFReadL(1, 10000, mem_f)
    vrt_template += data.decode('ascii')

  gdal.VSIFCloseL(mem_f)
  gdal.Unlink(vrt_tpl_file)

  return vrt_template

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="""Creating VRT to compute total_cover = nphot + phot on the fly""")
  parser.add_argument('src_ds', type=str, help='input dataset in the format of NETCDF:<nc file name>:sub_dataset')
  parser.add_argument('--dst_dir', default='./', type=str, help='VRT output file name')
  args = parser.parse_args()

  parts = args.src_ds.split(':')
  assert len(parts) == 3, 'input dataset must be in the format of NETCDF:<nc file name>:sub_dataset'
  src_filename = os.path.abspath(parts[1])

  vrt_template = get_vrt_template(args.src_ds)

  root = ET.fromstring(vrt_template)

  bands = root.findall('./VRTRasterBand')
  for ib, band_ele in enumerate(bands):
    tc_band = get_total_cover_band(ib+1, src_filename)
    tc_band_ele = ET.fromstring(tc_band)
    root.remove(band_ele)
    root.append(tc_band_ele)
    
  ele_tree = ET.ElementTree(root)

  base_src_filename = os.path.basename(src_filename)
  dst_filename = os.path.splitext(base_src_filename)[0] + '.vrt'
  dst_path = os.path.join(args.dst_dir, dst_filename)
  ele_tree.write(dst_path)
