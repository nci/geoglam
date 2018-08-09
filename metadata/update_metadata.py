import argparse
import json
import netCDF4

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Geoglam update metadata")
    parser.add_argument("--update_file", type=str, help="Existing prodcut file to update")
    parser.add_argument("--metadata_file", type=str, help="Metadata file")
    args = parser.parse_args()

    with netCDF4.Dataset(args.update_file, 'a', format='NETCDF4') as nc:
        with open(args.metadata_file) as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(nc, key, attrs[key])        
