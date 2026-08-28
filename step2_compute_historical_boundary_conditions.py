import numpy as np
import pandas as pd
from numpy import inf
import time
import datetime as dt
import xarray as xr
import math
from scipy.stats import skew
import pickle
import json
import sys
from shapely.geometry import Point, LineString
import rioxarray

from config import TRANSECT_DATA, WORKING_DICTS, ERA5_WAVES, dfN_i, WORKING_DIR, CONED_DTM
from step1_compile_historical_boundary_data import wl_file, wl_time_file, tide_file, tide_time_file


### load the saved WL and Tide data dictionaries ###

with open(WORKING_DICTS + wl_file, 'rb') as f:
    station_data = pickle.load(f)
with open(WORKING_DICTS + wl_time_file, 'rb') as f:
    station_time = pickle.load(f)
    
with open(WORKING_DICTS + tide_file, 'rb') as f:
    station_data_tides = pickle.load(f)
with open(WORKING_DICTS + tide_time_file, 'rb') as f:
    station_time_tides = pickle.load(f)
 
    
# load transects dataframe
df = pd.read_csv(TRANSECT_DATA)
    
# get list of unique wave stations from the transects
DFMids = df['DFMid'].unique()


# In[ ]:  

### Compute Wave Power and WL stats for all DFMid ###

# define funciton to convert ERA5 ordinal time steps to datetime64 
def OrdinalToDatetime(ordinal):
    plaindate = dt.date.fromordinal(int(ordinal))
    date_time = dt.datetime.combine(plaindate, dt.datetime.min.time())
    return np.datetime64(date_time + dt.timedelta(days=ordinal-int(ordinal)))

# Define vars
g = 9.81 #m/s2, gravity
d = 1022.7 #kg/m3, density (% sigma t is usually used at around 1.0227 but waves are on the scale of 1000kg/m3)
con = (d*(g**2))/(64*math.pi) #constant, kg*m2/m3*s4, for wave power calc

# Define uncertainties
hs_unc = 0.25
tm_unc = 1.3

# open waves file
waves = xr.open_dataset(ERA5_WAVES)
# get wave data timestamps
wave_t = list(waves.time.values)

# initialize empty lists for outputs of each DFMid station
CumWvPwrlist = []
WvPwrSkewlist = []
sum_err = []
TWL_range = []
TWL_skew = []
Tidal_range = []
NTR_95th = [] 
mask_indices = []
st = time.time()
wave_idx_used = []

current_county = None
current_file = None

for i in range(len(DFMids)):
    
    #county = "King" # replace with df.County[i] if dataframe contains sites from more than one county
    dfmid = DFMids[i]
    DFMlon = df[df['DFMid'] == DFMids[i]]['DFMlon'].reset_index(drop=True)[0]
    DFMlat = df[df['DFMid'] == DFMids[i]]['DFMlat'].reset_index(drop=True)[0]
    
    # get water levels and timestamps
    twls = np.concatenate(station_data[dfmid]) / 10000
    wltime = np.concatenate(station_time[dfmid]) 
    tides = np.concatenate(station_data_tides[dfmid]) / 10000
    tidetime = np.concatenate(station_time_tides[dfmid]) 

    # define cutoff date
    cutoff = np.datetime64('1941-01-01')
    # create boolean mask
    mask = wltime >= cutoff
    tmask = tidetime >= cutoff 
    # apply mask to arrays
    twls = twls[mask]
    tides = tides[tmask]
    wltime = wltime[mask]
    
    wlrange = np.nanpercentile(twls, 99) - np.nanpercentile(twls, 1)
    wlskew = skew(twls)
    
    tiderange = np.nanpercentile(tides, 99) - np.nanpercentile(tides, 1)
    
    TWL_range.append(np.round(wlrange, 2))
    TWL_skew.append(np.round(wlskew, 2))
    
    Tidal_range.append(np.round(tiderange,2))

    # get 95th percentile of positive NTRs
    NTR = twls - tides
    ntr_filt = [item for item in NTR if item > 0.0] # remove NTRs less than 0
    NTR_95th.append(np.nanpercentile(ntr_filt,95))
  
    percent = np.round(i / len(DFMids) * 100, 2)
    print(str(percent) + "% complete, " + f"Processing station: {dfmid}")
    
    # get the lat lon from wave file
    mat_lons = waves['lon'].values
    mat_lats = waves['lat'].values
     
    wltime = pd.to_datetime(wltime)

    # define water level threshold
    threshold = ((np.nanpercentile(tides, 99) + np.nanpercentile(tides, 1)) / 3) * 2
    # identify all waterlevels over threshold (meters)
    mask = twls > (threshold) # water level reaches toe threshold
    # Indices where condition is true
    indices = np.where(mask)[0]
    # Check whether previous value was also > threshold
    prev_mask = np.roll(mask, 1)
    prev_mask[0] = False  # first element has no previous value
    # 1 if previous was also > threshold, else 0 (only for matching indices)
    continuity = prev_mask[indices].astype(int)
    # get timestamps for the waterlevels > 0.5
    t_exceed = wltime[indices]
        
    ## now to get wave data station by station ##
    # Compute distance to all LUT stations
    dist = np.sqrt((mat_lons - DFMlon)**2 + (mat_lats - DFMlat)**2)
    # Find nearest station
    idx = np.argmin(dist)
    wave_idx_used.append(idx)
    # safety check
    if dist[idx] > 0.01:   # adjust tolerance (~0.01° ≈ 1 km)
        print(f"Warning: large distance match for {dfmid}: {dist[idx]}")
    hs = waves['Hs'][:,idx].values # get hs values
    tm = waves['Tm'][:,idx].values / np.timedelta64(1, 's') # get tm values and convert from timedelta to float

    ## filter waves for just times when water level is above threshold
     
    # indices in t_exceed where a new exceedance period starts
    starts = np.where(continuity == 0)[0]
    # corresponding ends
    ends = np.empty_like(starts)
    for j, s in enumerate(starts):
        if j < len(starts) - 1:
            ends[j] = starts[j + 1] - 1
        else:
            ends[j] = len(t_exceed) - 1
     # build (start_time, end_time) intervals
    intervals = [(t_exceed[s], t_exceed[e]) for s, e in zip(starts, ends)]
    intervals = np.array(intervals, dtype='datetime64[ns]')
    if len(intervals) > 0:
         
        starts = intervals[:, 0]
        ends   = intervals[:, 1]
        wave_t = np.asarray(wave_t)
        # index of rightmost start <= wave_t
        interval_idx = np.searchsorted(starts, wave_t, side='right') - 1
        inside = ((interval_idx >= 0) & (wave_t <= ends[interval_idx]))
        matching_indices = np.where(inside)[0]
        mask_indices.append(matching_indices)
         
        hs_masked = hs[matching_indices]
        tm_masked = tm[matching_indices]
         
        # calculate wave powers from filtered data
        OffWvPwrfilt = ((con*(hs_masked**2)*tm_masked)) # calculate wave powers (kW/m)
        cumulative = np.nansum(OffWvPwrfilt)
        CumWvPwrlist.append(cumulative)
         
        WvPwrErr = con * np.sqrt(4 * hs_masked**2 * tm_masked**2 * hs_unc**2 + hs_masked**4 * tm_unc**2 ) # propogate error through wave power equation
        ErrSum = np.nansum(WvPwrErr)
        sum_err.append(ErrSum)
         
    else:
        CumWvPwrlist.append(np.nan)
        sum_err.append(np.nan)

    OffWvPwr = ((con*(hs**2)*tm)) # calculate unfiltered wave power (kW/m)
    WvPwrSkew = np.round(skew(OffWvPwr, nan_policy='omit'), 3)
    WvPwrSkewlist.append(WvPwrSkew)

dfnew = pd.DataFrame(DFMids, columns=["DFMid"])

dfnew['WavePower_Cum'] = CumWvPwrlist
dfnew['WavePower_Skew'] = WvPwrSkewlist
dfnew['TWL_Range'] = TWL_range
dfnew['TWL_Skew'] = TWL_skew
dfnew['Tidal_Range'] = Tidal_range  
dfnew['NTR95'] = NTR_95th
dfnew['NTR95_per_TWLrange'] = np.array(NTR_95th) / np.array(TWL_range)

## add new computed variables to full transects dataframe, and save boundary conds dataframe checkpoint
df_merged = df.merge(dfnew, on='DFMid', how='left')
df_merged.to_csv(WORKING_DIR + '/' + dfN_i + '.csv', index=False)
df = df_merged

## calculate and save relative wavepower uncertainty
wp_error = np.round(np.nanmean(sum_err) / np.nanmean(CumWvPwrlist), 2)
parameters = {"relative_wave_power_uncertainty": float(wp_error)}
with open(WORKING_DIR + '/wp_uncertainty.json', "w") as f:
    json.dump(parameters, f, indent=4)

# In[ ]:

### calculate Bearing of the transects ###

def calculate_bearing(df):
    # Calculate differences in coordinates
    delta_easting = df['BluffToe_UTMX'] - df['BluffTop_UTMX']
    delta_northing = df['BluffToe_UTMY'] - df['BluffTop_UTMY']
    # Calculate bearing in radians (-pi to pi)
    bearing_rad = np.arctan2(delta_easting, delta_northing)
    # Convert to degrees (0-360)
    bearing_deg = np.round(np.degrees(bearing_rad) % 360 , 2)
    return bearing_deg
    
df['Bearing'] = calculate_bearing(df) # add bearing to data frame


# In[ ]:

### Calculate depths of closure ###

g = 9.81 #m/s2, gravity
meandep = [] # initialize list for closure depths
wave_idx_used = []

# get the lats/lons from wave file
mat_lons = waves['lon'].values
mat_lats = waves['lat'].values

st = time.time()

j = 0
for j in range(len(DFMids)): #for all wave file names within the table df
    st1 = time.time()
    DFMid = DFMids[j] # wave and water level site ID

    # get DFMid lat lons
    DFMlon = df.loc[df['DFMid'] == DFMid, 'DFMlon'].iloc[0]
    DFMlat = df.loc[df['DFMid'] == DFMid, 'DFMlat'].iloc[0]

    ## now to get wave data station by station ##
    # Compute distance to all LUT stations
    dist = np.sqrt((mat_lons - DFMlon)**2 + (mat_lats - DFMlat)**2)
    # Find nearest station
    idx = np.argmin(dist)
    wave_idx_used.append(idx)
    # safety check
    if dist[idx] > 0.01:   # adjust tolerance (~0.01° ≈ 1 km)
        print(f"Warning: large distance match for {dfmid}: {dist[idx]}")
    hs = waves['Hs'][:,idx].values # get the wave heights at that index
    tm = waves['Tm'][:,idx].values / np.timedelta64(1, 's')  # get the wave periods at that index  
    valid = (tm > 0) & (hs > 0) & np.isfinite(tm) & np.isfinite(hs) # mask for only valid hs and tm values
    dep = 2.28 * hs[valid] - (68.5 * (hs[valid]**2/(g*tm[valid]**2))) # depth of closure, from Limber et al 
    depm = np.nanmean(dep) # mean depth of closure 
    meandep.append(depm) # append mean closure depth to list

    et = time.time()
    sys.stdout.write('\rElapsed time: %.ds; Sites Completed: %.d; Remaining Sites: %s' %((et-st),(j+1),(len(DFMids)-(j+1))))
    sys.stdout.flush()

dfdep = pd.DataFrame(columns = ['DFMid', 'ClosureDepth'])
dfdep['DFMid'] = DFMids
dfdep['DepthClosure_meters'] = pd.to_numeric(meandep)


## add computed closure depths to the full boundary conditions dataframe, and save boundary conds dataframe checkpoint

df_merged2 = df.merge(dfdep, on='DFMid', how='left')
df_merged2.to_csv(WORKING_DIR + '/' + dfN_i + '.csv', index=False)
df = df_merged2


# In[ ]:

### find the point along transect where DTM elevation == depth of closure ###

def print_progress(iteration, total):
    current = "{:.0f}".format(iteration)
    print(f'\rProcessing: {current} of '+str(total)+'\n', end="", flush=True)

# NOTE!!! DTM and the bluff toe points must be in same CRS and lat/lon formatted as easting/northing in meters (they should be if using the published data)
 
def create_line_from_bearing( lon, lat, bearing, distance=500):
    """
    Creates a line starting from `point` in the direction of `bearing` for a specified `distance` (in meters).
    """
    # Convert bearing to radians
    bearing_rad = math.radians(bearing)
    
    # create point
    start = Point(lon, lat)
    
    # Calculate the end point
    end_x = lon + distance * math.sin(bearing_rad)
    end_y = lat + distance * math.cos(bearing_rad)
    
    # Create the line
    return LineString([start, Point(end_x, end_y)])

lines_list = []
# Add the lines to the GeoDataFrame
lines_list = df.apply(lambda row: create_line_from_bearing(row.BluffToe_UTMX, row.BluffToe_UTMY, row.Bearing), axis=1)

# Load xarray
dtm = rioxarray.open_rasterio(CONED_DTM)  # open geotiff raster again as xarray

## find 'closure depth' point, seaward of toe ##   
def find_closure_depth_point(line, dep, dtm, step=1):
    """
    Samples along the line to find the first point where the DTM value is less than or equal to closure depth
    """
    
    # Get the length of the line in meters
    line_length = line.length
    
    # Sample points along the line at the given step size
    num_samples = int(line_length / step)
    sampled_points = [line.interpolate(step * i) for i in range(num_samples + 1)]
    
    for point in sampled_points:
        # Get the DTM value at the point
        value = dtm.sel(x=point.x, y=point.y, method="nearest").values
        
        # Check if the value is less than or equal to closure depth
        if value <= dep:
            return point
        
    return None  # Return None if no point meets the condition

# Apply the function to each row to find the closure depth point
dpoint_list = []
for i in range(0, len(lines_list)):
    dep = df['DepthClosure_meters'][i]
    dpoint_list.append(find_closure_depth_point(lines_list[i], dep, dtm))
    et = time.time()
    sys.stdout.write('\rElapsed time: %.ds; Sites Completed: %.d; Remaining Sites: %s' %((et-st),(i+1),(len(lines_list)-(i+1))))
    sys.stdout.flush()


# In[ ]:
    
## find active beach width ##
def calculate_distance(cdepth, toe):
    """
    Calculates the distance between the original point (toe_lon, toe_lat) and the closure depth point.
    """
    # Check if the cdepth_point exists
    if cdepth is None:
        return 50
    
    # Calculate the distance (in meters, since CRS is in meters)
    distance = toe.distance(cdepth)
    
    return distance

# Apply the distance calculation function
width = []
for i in range(0, len(dpoint_list)):
    toe = Point(df.BluffToe_UTMX[i], df.BluffToe_UTMY[i])
    w = calculate_distance(dpoint_list[i], toe)
    if w > 50:
        width.append(50)
    else:
        width.append(w)
 
df['BeachWidth_meters'] = width


# In[ ]:

## find beach slope ##
rise = df.BluffToe_Elev - df.DepthClosure_meters
run = pd.Series(width)
beach_slope_pct = rise/run * 100

df['ActiveSlope_pct'] = beach_slope_pct


# In[ ]:
    
### save the new boundary conditions file with all the computed values to a csv file ###

df = df.drop_duplicates(subset=["OutputID"])
df.to_csv(WORKING_DIR + '/' + dfN_i + '.csv', index=False)