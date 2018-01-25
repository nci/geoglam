This script downloads and transforms CHIRPS2.0 from the FTP server to the local destination on disk.

To run this script:

1.- This script requires specific versions of GDAL and NetCDF4 which are available on the NCI Raijin and VDI nodes.

2.- Make sure the script is executable:
	chmod +x chirps.sh

3.- Pass the year as an argument to the script:
	./chirps 2017

This script can be rerun many times for the same year as a way of updating the contents of the file when new data is available for the current year.

National Computational Infrastructure (ANU) 2017
