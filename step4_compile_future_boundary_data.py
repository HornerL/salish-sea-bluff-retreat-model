import pandas as pd
import numpy as np
import xarray as xr
import os
import time
import pickle
from collections import defaultdict

from config import WORKING_DIR, WORKING_DICTS, TRANSECT_DATA, hist_rates_file, SLR, CMIP6WAVES, CMIP6WL_DIR, CMIP6TIDE_DIR, CMIPDIFF_DIR

### load historic bootstrap predictions, if you want to skip the previous steps
df_mvar = pd.read_csv(WORKING_DIR + '/' + hist_rates_file + '.csv')

### fix potential issue with column data type 
df_mvar = df_mvar.replace('#NAME?', np.nan)
df_mvar['ActiveSlope_pct'] = pd.to_numeric( df_mvar['ActiveSlope_pct'])

### convert slope from percent to unitless ###
beach_slope = (df_mvar['ActiveSlope_pct']/100)
df_mvar['Beach_Slope'] = beach_slope

### define additional variable measurement errors
tm_err = 1.3
hs_err = 0.25
wl_err = 0.15
elev_err = 0.25
slope_err = ((elev_err * np.sqrt(2)) / df_mvar['BeachWidth_meters']) # calculate uncertainty of slope


# In[ ]:
    
### Now we get the boundary conditions from future CMIP6 data

# Load File containing all transects
filename = TRANSECT_DATA
df2 = pd.read_csv(filename)

# get list of unique wave stations from the transects
DFMids = df2['DFMid'].unique()
stations_bytes = {s: s.encode("utf-8") for s in DFMids}

## CMIP6 waves file
filepath_wv6 = CMIP6WAVES + SLR +"m.nc"

## directory containing CMIP6 wl files
directory_wl6 = CMIP6WL_DIR + SLR 

## directory containing SLRxxx tides files
directory_td6 = CMIP6TIDE_DIR + '/' + SLR

## directory containing cmip_diff files
directory_cdiff = CMIPDIFF_DIR + '/' + SLR

## Get CMIP6 waves
waves6 = xr.open_dataset(filepath_wv6)

## Get a list of all CMIP6 wl files
filepathswl6 = [os.path.join(directory_wl6, f) for f in os.listdir(directory_wl6) if os.path.isfile(os.path.join(directory_wl6, f))]

## Get a list of all CMIP6 tide files
filepathstd6 = [os.path.join(directory_td6, f) for f in os.listdir(directory_td6) if os.path.isfile(os.path.join(directory_td6, f))]

## get list of all CMIP_Diff wl files
filepathscdiff = [os.path.join(directory_cdiff, f) for f in os.listdir(directory_cdiff) if os.path.isfile(os.path.join(directory_cdiff, f))]



## loop over the cmipdiff wl files, average the cmipdiffs, and build an array for each DFMid ##

cmipavg = defaultdict(list)
cmip_time = defaultdict(list)

DFMids_list = list(DFMids)

cmipavg = defaultdict(list)
cmip_time = defaultdict(list)


print("Looping through cmipdiff files, loading data\n")
for i, pathcdiff in enumerate(filepathscdiff):

    et = time.time()
    print(f"\rProcessing file {i+1}/{len(filepathscdiff)}", end="", flush=True)

    with xr.open_dataset(pathcdiff, chunks={'time': 500}) as cmip:

        wl_stations = cmip.coords['station'].values
        dtime = cmip['time'].values
    
        station_to_idx = {s: j for j, s in enumerate(wl_stations)}
        valid = [(s, station_to_idx[s]) for s in DFMids_list if s in station_to_idx] ## makes sure the station is in the list of King stations
    
        if not valid:
            continue
    
        stations, indices = zip(*valid)
    
        subset = cmip['cmip_diff'].isel(station=list(indices))
    
        ## mask bad values
        fill_val = cmip['cmip_diff'].attrs.get('_FillValue', -2147483647)
        subset = subset.where(subset != fill_val)
    
        ## mean across cmip6
        mean_vals = subset.mean(dim='cmip6').load()
    
        for k, station in enumerate(stations):
            cmipavg[station].append(mean_vals[:, k].values / 10000)
            cmip_time[station].append(dtime)

### save cmipdiff data 

with open(WORKING_DICTS + "/Cmipdiff_wl_SLR" + (SLR) + "_defaultdict.pkl", 'wb') as f:
    pickle.dump(cmipavg, f)


# In[ ]:
## loop over the water level files and build an array for each DFMid ##

wl_data = defaultdict(list)
wl_time = defaultdict(list)

DFMids_list = list(DFMids)

print(" Looping through WL files, loading data\n")

for i, pathwl in enumerate(filepathswl6):

    print(f"\rProcessing file {i+1}/{len(filepathswl6)}", end="", flush=True)

    try:
        with xr.open_dataset(pathwl, chunks={'time': 1000}) as ds:

            # Filter out data before 1941
            #ds = ds.sel(time=slice("1941-01-01", None))
            
            wl_stations = [s.decode('utf-8') if isinstance(s, bytes) else s for s in ds['station'].values]
            dtime = ds['time'].values

            # build lookup
            station_to_idx = {s: j for j, s in enumerate(wl_stations)}

            valid = [(s, station_to_idx[s]) for s in DFMids_list if s in station_to_idx]

            if not valid:
                continue

            stations, indices = zip(*valid)

            # subset only needed stations
            subset = ds['waterlevel'].isel(station=list(indices))

            # load into memory AFTER subsetting
            subset = subset.load()   # (time, selected_stations)

            # store per station
            for k, station in enumerate(stations):
                wl_data[station].append(subset[:, k].values)
                wl_time[station].append(dtime)

    except Exception as e:
        print(f"\nSkipping file {i}: {e}")
        continue

with open(WORKING_DICTS + "/ERA5_WL_SLR" + (SLR) + "_defaultdict.pkl", 'wb') as f:
    pickle.dump(wl_data, f)
with open(WORKING_DICTS + "/ERA5_WLtime_SLR" + (SLR) + "_defaultdict.pkl", 'wb') as f:
    pickle.dump(wl_time, f)
    

# In[ ]:    
## loop over the tides files and build an array for each DFMid ##

tide_data = defaultdict(list)
tide_time = defaultdict(list)

DFMids_list = list(DFMids)

print(" Looping through WL files, loading data\n")

for i, pathtd in enumerate(filepathstd6):

    print(f"\rProcessing file {i+1}/{len(filepathstd6)}", end="", flush=True)

    try:
        with xr.open_dataset(pathtd, chunks={'time': 1000}) as ds:

            ## Filter out data before 1941
            #ds = ds.sel(time=slice("1941-01-01", None))

            td_stations = [s.decode('utf-8') if isinstance(s, bytes) else s for s in ds['station'].values]
            dtime = ds['time'].values

            ## build lookup
            station_to_idx = {s: j for j, s in enumerate(td_stations)}

            valid = [(s, station_to_idx[s]) for s in DFMids_list if s in station_to_idx]

            if not valid:
                continue

            stations, indices = zip(*valid)

            ## subset only needed stations
            subset = ds['waterlevel'].isel(station=list(indices))

            ## load into memory AFTER subsetting
            subset = subset.load()   ## (time, selected_stations)

            ## store per station
            for k, station in enumerate(stations):
                tide_data[station].append(subset[:, k].values)
                tide_time[station].append(dtime)

    except Exception as e:
        print(f"\nSkipping file {i}: {e}")
        continue

with open(WORKING_DICTS + "/ERA5_tidetime_SLR" + (SLR) + "_defaultdict.pkl", 'wb') as f:
    pickle.dump(tide_time, f)
with open(WORKING_DICTS + "/ERA5_tide_SLR" + (SLR) + "_defaultdict.pkl", 'wb') as f:
    pickle.dump(tide_data, f)
    

# In[ ]:    
## loop over the water level files and build a lat/lon array for each DFMid ##

wl_lon = defaultdict(list)
wl_lat = defaultdict(list)

print("\n Looping through TWL files, loading data\n")
for j, filepathwl in enumerate(filepathswl6):
    print(f"\rProcessing file {j+1}/{len(filepathswl6)}", end="", flush=True)
    with xr.open_dataset(filepathwl) as wl:
        wl_stations = wl['station'].values

        for dfmid, station in stations_bytes.items():
            if station in wl_stations:
                da = wl.sel(station=station)
                wl_lon[dfmid].append(da['lon'].values)
                wl_lat[dfmid].append(da['lat'].values)

with open(WORKING_DICTS + "/ERA5_WLlat_SLR" + (SLR) + "_defaultdict.pkl", 'wb') as f:
    pickle.dump(wl_lat, f)
with open(WORKING_DICTS + "/ERA5_WLlon_SLR" + (SLR) + "_defaultdict.pkl", 'wb') as f:
    pickle.dump(wl_lon, f)
    
    
