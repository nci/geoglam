# Geoglam Fractional cover - MODIS, CSIRO Land and Water algorithm

Vegetation fractional cover represents the exposed proportion of Photosynthetic Vegetation (PV), Non-Photosynthetic Vegetation (NPV) and Bare Soil (BS) within each pixel. In forested canopies the photosynthetic or non-photosynthetic portions of trees may obscure those of the grass layer and/or bare soil. The MODIS Fractional Cover product is derived from the MODIS Nadir BRDF-Adjusted Reflectance (NBAR) product (MCD43A4, collection 5). A suite of derivative are also produced, namely total vegetation cover (PV+NPV), monthly fractional cover and total vegetation cover, monthly anomaly of total cover against the time series, and three-monthly total cover difference. MODIS fractional cover has been validated for Australia.

### Enivronment setup

* Run `git clone https://github.com/nci/geoglam.git` into your local directory
* Run `setup_env.sh` (i.e. `bash setup_env.sh`) to set up miniconda environment. You only need to run this step once
* Open `submit_pbs_jobs.sh` and fill up `NCI_PROJECT` and `BASE_DIR`. `NCI_PROJECT` is your NCI project number, `fr1`, for example. `BASE_DIR` is the full path of the to-be-generated geoglam products. For example, `/short/<project no.>/<user id>/geoglam_data`

### Generating Geoglam netCDF products 

* Once the miniconda environment is set up, run `submit_pbs_jobs.sh <tile csv file>` to submit PBS jobs for each tile in the tile csv file. There are a few example tile csv files under the `tile_files` directory.
* Once the PBS jobs finished, `run utils/parse_logs.py` (i.e. `python utils/parse_logs.py <log directory>`) to extract any failed tiles. If there are failed tiles, a csv file that contains the failed tiles will be written under the current working directory.
* You might want to archive the log files for future references. To do so, run `utils/archive_log_files.sh`

### Updating netCDF metadata

If the metadata of the netCDF product files need to be updated, please do the following steps:
* Enter the `metadata` directory (i.e. `cd metadata`) 
* Run `bash batch_update_metadata.sh <metadata json file>`. The `<metadata json file>` has the same format as `nc_metadata.json` located under the `geoglam` root directory.

Optionally, one might want to open `batch_update_metadata.sh` and edit the `UPDATE_DIR` and the `PYTHON` variables to point to alternative locations.
