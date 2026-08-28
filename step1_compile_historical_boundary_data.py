### Only run this code one time per set of transects ###

import pandas as pd
from collections import defaultdict
import numpy as np
import xarray as xr
import os
import pickle

from config import TRANSECT_DATA, ERA5_WL_DIR, ERA5_TIDE_DIR, WORKING_DICTS

# filnames for time series working dicts
wl_file = "/ERA5_WL_SLR000_defaultdict.pkl"
wl_time_file = "/ERA5_WLtime_SLR000_defaultdict.pkl" 
tide_file = "/ERA5_tide_SLR000_defaultdict.pkl"
tide_time_file = "/ERA5_tidetime_SLR000_defaultdict.pkl"

def compile_data():
        
    # load transects dataframe
    df = pd.read_csv(TRANSECT_DATA)
    
    # get list of unique wave stations from the transects
    DFMids = df['DFMid'].unique()
    stations_bytes = {s: s.encode("utf-8") for s in DFMids}
    
    
    ### Load the WL and tides data and build an array for each DFMID, then save the compiled dictionaries ###
    
    # Get a list of all wl files
    filepaths1 = [os.path.join(ERA5_WL_DIR, f) for 
                  f in os.listdir(ERA5_WL_DIR) if os.path.isfile(os.path.join(ERA5_WL_DIR, f))]
    
    ## loop over the water level files and build an array for each DFMid ##
    
    station_data = defaultdict(list)
    station_time = defaultdict(list)
    
    print("Looping through TWL files, loading data\n")
    
    for i, filepath in enumerate(filepaths1):
    
        print(f"\rProcessing file {i+1}/{len(filepaths1)}", end="", flush=True)
    
        try:
            with xr.open_dataset(filepath, chunks={'time': 1000}) as wl:
    
                wl_stations = wl['station'].values
                dtime = wl['time'].values
    
                # Build lookup: station name/bytes -> station index
                station_to_idx = {station: j for j, station in enumerate(wl_stations)}
    
                # Find only the DFMids/stations that are present in this file
                valid = [(dfmid, station_to_idx[station])
                    for dfmid, station in stations_bytes.items()
                    if station in station_to_idx]
    
                if not valid:
                    continue
    
                dfmids, indices = zip(*valid)
    
                # Subset only the needed stations before loading into memory
                subset = wl['waterlevel'].isel(station=list(indices))
    
                # Load the reduced array once instead of doing .sel() repeatedly
                subset = subset.load()
    
                # Store each station
                for k, dfmid in enumerate(dfmids):
    
                    station_data[dfmid].append(subset[:, k].values.astype(np.float32))
    
                    station_time[dfmid].append(dtime)
    
        except Exception as e:
            print(f"\nSkipping file {i+1}: {e}")
            continue
    
    
    ## loop over the tides files and build an array for each DFMid ##
    
    # Get a list of all tide wl files
    filepaths2 = [os.path.join(ERA5_TIDE_DIR, f) for f in os.listdir(
        ERA5_TIDE_DIR) if os.path.isfile(os.path.join(ERA5_TIDE_DIR, f))]
    
    station_data_tides = defaultdict(list)
    station_time_tides = defaultdict(list)
    
    print("\nLooping through Tides files, loading data\n")
    
    for i, filepath in enumerate(filepaths2):
    
        print(f"\rProcessing file {i+1}/{len(filepaths2)}", end="", flush=True)
    
        try:
            with xr.open_dataset(filepath, chunks={'time': 1000}) as wl:
    
                wl_stations = wl['station'].values
                dtime = wl['time'].values
    
                # Build lookup: station -> station index
                station_to_idx = {station: j for j, station in enumerate(wl_stations)}
    
                # Find only the DFMids/stations present in this file
                valid = [(dfmid, station_to_idx[station])
                    for dfmid, station in stations_bytes.items()
                    if station in station_to_idx]
    
                if not valid:
                    continue
    
                dfmids, indices = zip(*valid)
    
                # Subset only the needed stations
                subset = wl['waterlevel'].isel(station=list(indices))
    
                # Load into memory only AFTER subsetting
                subset = subset.load()
    
                # Store per DFMid
                for k, dfmid in enumerate(dfmids):
    
                    station_data_tides[dfmid].append(subset[:, k].values.astype(np.float32))
    
                    station_time_tides[dfmid].append(dtime)
    
        except Exception as e:
            print(f"\nSkipping file {i+1}: {e}")
            continue
        
    # Save wl and tide dictionaries 
    with open(WORKING_DICTS + wl_file, 'wb') as f:
        pickle.dump(station_data, f)
    with open(WORKING_DICTS + wl_time_file, 'wb') as f:
        pickle.dump(station_time, f)
    with open(WORKING_DICTS + tide_file, 'wb') as f:
        pickle.dump(station_data_tides, f)
    with open(WORKING_DICTS + tide_time_file, 'wb') as f:
        pickle.dump(station_time_tides, f)
        
if __name__ == "__main__":
    compile_data()
