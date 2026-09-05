import pickle
import pandas as pd
import numpy as np
import xarray as xr
import psutil
import os
import math
import time
from scipy.stats import skew

from config import WORKING_DIR, WORKING_DICTS, SLR, CMIP6WAVES, dfN_i, TRANSECT_DATA

## get list of unique wl stations from the transects
df = pd.read_csv(TRANSECT_DATA)
DFMids = df['DFMid'].unique()

### Load waves, WL and Tide data here if you have already compiled them for the given SLR

## CMIP6 waves file
filepath_wv6 = CMIP6WAVES + SLR +"m.nc"

## process data in chunks to work around potential RAM limitations with large counties
n_chunks = 4
for k in range (0,4):
    
    with open(WORKING_DICTS + "/ERA5_WL_SLR" + SLR + "_defaultdict.pkl", 'rb') as f:
        wl_data = pickle.load(f)
       
    #all_stations = list(wl_data.keys())
    station_chunks = np.array_split(DFMids, n_chunks)   
    chunk = station_chunks[k] 
        
    wl_data = {k: wl_data[k] for k in chunk}
            
    with open(WORKING_DICTS + "/ERA5_WLtime_SLR" + SLR + "_defaultdict.pkl", 'rb') as f:
        wl_time = pickle.load(f)
        
    wl_time = wl_time[chunk[0]]
    
    with open(WORKING_DICTS + "/ERA5_tide_SLR" + SLR + "_defaultdict.pkl", 'rb') as f:
        tide_data = pickle.load(f)
        
    tide_data = {k: tide_data[k] for k in chunk}
    
    with open(WORKING_DICTS + "/ERA5_tidetime_SLR" + SLR + "_defaultdict.pkl", 'rb') as f:
        tide_time = pickle.load(f)
        
    tide_time = tide_time[chunk[0]]
        
    with open(WORKING_DICTS + "/ERA5_WLlon_SLR" + SLR + "_defaultdict.pkl", 'rb') as f:
        wl_lon = pickle.load(f)
    
    for dfmid in wl_lon:
        wl_lon[dfmid] = wl_lon[dfmid][0]
       
    
    with open(WORKING_DICTS + "/ERA5_WLlat_SLR" + SLR + "_defaultdict.pkl", 'rb') as f:
        wl_lat = pickle.load(f)
    
    for dfmid in wl_lat:
        wl_lat[dfmid] = wl_lat[dfmid][0]
        
       
    with open(WORKING_DICTS + "/Cmipdiff_wl_SLR" + SLR + "_defaultdict.pkl", 'rb') as f:
        cmipavg = pickle.load(f)
    
    cmipavg = {k: cmipavg[k] for k in chunk}
    
    
    with open(WORKING_DICTS + "/ERA5_tide_SLR000_defaultdict.pkl", 'rb') as f:
        station_data_tides = pickle.load(f)
    with open(WORKING_DICTS + "/ERA5_tidetime_SLR000_defaultdict.pkl", 'rb') as f:
        station_time_tides = pickle.load(f)
    
    station_data_tides = {k: station_data_tides[k] for k in chunk}
    
    station_time_tides = station_time_tides[chunk[0]]
    
    
    ## find needed wave stations
    
    waves6 = xr.open_dataset(filepath_wv6)
    
    wave_lon = waves6['lon'].values
    wave_lat = waves6['lat'].values
    
    needed_wave_idx = set()
    
    for dfmid in chunk:
    
        DFMlon = wl_lon[dfmid].item()
        DFMlat = wl_lat[dfmid].item()
    
        dist = np.sqrt((wave_lon - DFMlon)**2 +
                       (wave_lat - DFMlat)**2)
    
        idx = np.argmin(dist)
    
        needed_wave_idx.add(idx)
    
    needed_wave_idx = np.array(sorted(needed_wave_idx))
    
    waves6 = waves6.isel(station=needed_wave_idx)
    
    ## filename that the computed SLRxxx boundary conditions will be saved to
    dfN_SLR = dfN_i + '_SLR' + SLR
    
    ### compute wave and wl boundary conditions for SLRxxx scenario ###
    
    ## Define vars
    g = 9.81 ## m/s2, gravity
    d = 1022.7 ## kg/m3, density (% sigma t is usually used at around 1.0227 but waves are on the scale of 1000kg/m3)
    con = (d*(g**2))/(64*math.pi) ## constant, kg*m2/m3*s4, for wave power calc
    
    ## Define uncertainties
    hs_unc = 0.25
    tm_unc = 1.3
    
    ## initialize empty lists for outputs of each DFMid station
    CumWvPwrlist = []
    WvPwrSkewlist = []
    TWL_range = []
    TWL_skew = []
    NTR_95th = [] 
    sum_err = []
    Tidal_range = []
    mask_indices = []
    wave_idx_used = []
    st = time.time()
    
    wave_lon = waves6['lon'].values
    wave_lat = waves6['lat'].values
    wave_dtime = waves6['time'].values
      
    for i in range(len(chunk)):
        
        import gc
        gc.collect()
    
        if i % 50 == 0:
            process = psutil.Process(os.getpid())
            print(f"Memory = {process.memory_info().rss/1024**3:.2f} GB" )
        
        dfmid = chunk[i]
    
        slr = (pd.to_numeric(SLR)) / 100  ## convert slr to meters
        
        ## get water levels and timestamps
        wls = np.concatenate(wl_data[dfmid]) / 10000 
        cdiffs = np.concatenate(cmipavg[dfmid]) 
        wltime = np.concatenate(wl_time) 
        tides = np.concatenate(station_data_tides[dfmid]) / 10000 ## historic tides
        tidetime = np.concatenate(station_time_tides) 
        tides2 = np.concatenate(tide_data[dfmid]) / 10000 ## future tides
        #cmtime = np.concatenate(cmip_time[dfmid])
        tidetime2 = np.concatenate(tide_time)
        twls = wls + cdiffs - slr ## Calculate the new wl values 
        
        
        ## define cutoff date
        cutoff = np.datetime64('1945-01-01')  ########## change date to 1945 for slr300 ############
        ## create boolean mask
        mask = wltime >= cutoff
        tmask = tidetime >= cutoff
        tmask2 = tidetime2 >= cutoff 
        ## apply mask to arrays
        twls = twls[mask]
        tides = tides[tmask]
        tides2 = tides2[tmask2]
        wltime = wltime[mask]
        
        wlrange = np.nanpercentile(twls, 99) - np.nanpercentile(twls, 1)
        wlskew = skew(twls, nan_policy='omit')
        
        tiderange = np.nanpercentile(tides2, 99) - np.nanpercentile(tides2, 1)
        
        TWL_range.append(np.round(wlrange, 2))
        TWL_skew.append(np.round(wlskew, 2))
        
        Tidal_range.append(np.round(tiderange,2))
    
        ## get 95th percentile of positive NTRs
        NTR = twls - tides2
        ntr_filt = [item for item in NTR if item > 0.0] ## remove NTRs less than 0
        NTR_95th.append(np.nanpercentile(ntr_filt,95))
            
        wltime = pd.to_datetime(wltime)
    
        ## define water level threshold
        threshold = ((np.nanpercentile(tides, 99) + np.nanpercentile(tides, 1)) / 3) * 2
        ## identify all waterlevels over threshold (meters)
        mask = twls > (threshold) ## water level reaches toe threshold
        ## Indices where condition is true
        indices = np.where(mask)[0]
        ## Check whether previous value was also > threshold
        prev_mask = np.roll(mask, 1)
        prev_mask[0] = False  ## first element has no previous value
        ## 1 if previous was also > threshold, else 0 (only for matching indices)
        continuity = prev_mask[indices].astype(int)
        ## get timestamps for the waterlevels > 0.5
        t_exceed = wltime[indices]
            
        ## now to get wave data station by station ##
        DFMlon = wl_lon[dfmid].item()
        DFMlat = wl_lat[dfmid].item()
        
        ## Compute distance to all LUT stations
        dist = np.sqrt((wave_lon - DFMlon)**2 + (wave_lat - DFMlat)**2)
        ## Find nearest station
        idx = np.argmin(dist)
        #wave_idx_used.append(idx)
        ## safety check
        if dist[idx] > 0.01:   ## adjust tolerance (~0.01° ≈ 1 km)
            print(f"Warning: large distance match for {dfmid}: {dist[idx]}")
        hs_cmipavg = waves6['hs_CmipDiff'].mean(dim='cmip6')[:,idx].values # average the hs cmip6 diff values
        hs = waves6['Hs'][:,idx].values + hs_cmipavg ## get the wave heights at that index
        tm = waves6['Tm'][:,idx].values #/ np.timedelta64(1, 's') # get tm values and convert from timedelta to float
    
        ## filter waves for just times when water level is above threshold
         
        ## indices in t_exceed where a new exceedance period starts
        starts = np.where(continuity == 0)[0]
        ## corresponding ends
        ends = np.empty_like(starts)
        for j, s in enumerate(starts):
            if j < len(starts) - 1:
                ends[j] = starts[j + 1] - 1
            else:
                ends[j] = len(t_exceed) - 1
        ## build (start_time, end_time) intervals
        intervals = [(t_exceed[s], t_exceed[e]) for s, e in zip(starts, ends)]
        intervals = np.array(intervals, dtype='datetime64[ns]')
        if len(intervals) > 0:
             
            starts = intervals[:, 0]
            ends   = intervals[:, 1]
            ## index of rightmost start <= wave_t
            interval_idx = np.searchsorted(starts, wave_dtime, side='right') - 1
            inside = ((interval_idx >= 0) & (wave_dtime <= ends[interval_idx]))
            matching_indices = np.where(inside)[0]
            #mask_indices.append(matching_indices)
             
            hs_masked = hs[matching_indices]
            tm_masked = tm[matching_indices]
             
            ## calculate wave powers from filtered data
            OffWvPwrfilt = ((con*(hs_masked**2)*tm_masked)) ## calculate wave powers (kW/m)
            cumulative = np.nansum(OffWvPwrfilt)
            CumWvPwrlist.append(np.round(cumulative,2))
             
            WvPwrErr = con * np.sqrt(4 * hs_masked**2 * tm_masked**2 * hs_unc**2 + hs_masked**4 * tm_unc**2 ) # propogate error through wave power equation
            ErrSum = np.round(np.nansum(WvPwrErr),1)
            sum_err.append(ErrSum)
             
        else:
            CumWvPwrlist.append(np.nan)
            sum_err.append(np.nan)
    
        OffWvPwr = ((con*(hs**2)*tm)) ## calculate unfiltered wave power (kW/m)
        WvPwrSkew = np.round(skew(OffWvPwr, nan_policy='omit'), 3)
        WvPwrSkewlist.append(WvPwrSkew)
        
        percent = np.round(i / len(chunk) * 100, 2)
        print(str(percent) + "% complete, " + f"Processing station: {dfmid}")
    
    dfnew2 = pd.DataFrame(chunk, columns=["DFMid"])
    
    dfnew2['WavePower_Cum'] = CumWvPwrlist
    dfnew2['WavePower_Skew'] = WvPwrSkewlist
    dfnew2['WavePower_Cum_err'] = sum_err 
    dfnew2['TWL_Range'] = TWL_range
    dfnew2['TWL_Skew'] = TWL_skew
    dfnew2['Tidal_Range'] = Tidal_range  
    dfnew2['NTR95'] = NTR_95th
    dfnew2['NTR95_per_TWLRange'] = np.array(NTR_95th) / np.array(TWL_range)
    
    
    ## modify NTR for SLR 
    
    slr = int(SLR) / 100
    dfnew2['NTR95'] = dfnew2['NTR95'] + slr
    dfnew2['NTR95_per_TWLRange'] = (dfnew2['NTR95']) / dfnew2['TWL_Range']
    
    dfnew2.to_csv(WORKING_DIR + '/CMIP6_SLR'+ SLR +'_WP_Chunk'+ str(k) +'.csv')


# In[ ]:

# combine chunked df into single df 

df_chunk3 = pd.read_csv(WORKING_DIR + '/CMIP6_SLR'+ SLR +'_WP_Chunk3.csv')
df_chunk2 = pd.read_csv(WORKING_DIR + '/CMIP6_SLR'+ SLR +'_WP_Chunk2.csv')
df_chunk1 = pd.read_csv(WORKING_DIR + '/CMIP6_SLR'+ SLR +'_WP_Chunk1.csv')
df_chunk0 = pd.read_csv(WORKING_DIR + '/CMIP6_SLR'+ SLR +'_WP_Chunk0.csv')

dfnew2 = pd.concat([df_chunk0, df_chunk1, df_chunk2, df_chunk3], ignore_index=True) 

## load a copy of the original boundary conditions and the computed historic boundary conditions
df2 = pd.read_csv(TRANSECT_DATA)
df_mvar = pd.read_csv(WORKING_DIR + '/' + dfN_i + '.csv')

## add new computed variables to transects dataframe

df_merged = df2.merge(dfnew2, on='DFMid', how='left')
df_merged.to_csv(WORKING_DIR + '/' + dfN_SLR + '.csv', index=False)

df2 = df_merged.copy()

df_mvar = df_mvar.sort_values(by='OutputID')
df2 = df2.sort_values(by='OutputID')
df2['Bearing'] = df_mvar['Bearing']

###### This saved file will contain ALL boundary conditions for SLRXXX Co Transects ######

df2 = df2.drop_duplicates(subset=["OutputID"])
df2.to_csv(WORKING_DIR + '/' + dfN_SLR + '.csv', index=False)
