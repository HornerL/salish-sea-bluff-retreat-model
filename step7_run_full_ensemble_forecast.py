# ## Now we apply TWL model, Hackney model, and Bruun Rule to estimate SLR impact (then average with the MV regression results)
# 
# For a given SLR secnario, we estimate the new erosion rate at ~20 yr intervals until the target date
# 
# Use known measurment error and basic error propogation rules for each model to find upper and lower uncertainty
# 
# Calculate weighted average of the FOUR models (Bruun, Hackney, TWL, and our future MV regression) using Limber's method (model_weight = 1/(Rd + 2*model_uncertainty)^2. Where Rd is the difference between the modeled rate and the mean of the other three models)
# 
# We run the models three times, once with our best estimate EPR rate, once with our lower prediction interval, and once with out upper prediction interval. Add the model uncertainty to the upper prediction interval rate, and subtract it from the lower prediction interval rate. 

import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import pickle
import math
import time
import sys
import os

from config import WORKING_DIR, WORKING_DICTS, SLR, ERA5_WAVES, TRANSECT_DATA, county, startyear, stopyear, hist_rates_file

dfN_Fcast = 'AllForecastRates_' + county + '_' + SLR + 'cm_by' + stopyear
dfN_Fcast2 = 'JustEnsembleRates_' + county + '_' + SLR + 'cm_by' + stopyear

df2 = pd.read_csv(TRANSECT_DATA)
df_mvar = pd.read_csv(WORKING_DIR + '/' + hist_rates_file + '.csv')
df_mvar2 = pd.read_csv(WORKING_DIR + '/CMIP6_SLR' + SLR +'_MVrates.csv')

## make sure all dataframes are in the same order

df_mvar = df_mvar.sort_values(by='OutputID')
df_mvar2 = df_mvar2.sort_values(by='OutputID')
df2 = df2.sort_values(by='OutputID')

## insure column header consistency
df_mvar['BluffToe_UTMY'] = df2['BluffToe_UTMY'] 
df_mvar['BluffToe_UTMX'] = df2['BluffToe_UTMX'] 
df_mvar['BluffToe_Elev'] = df2['BluffToe_Elev'] 


## load the saved WL and Tide data for SLR000
# filnames for time series working dicts
wl_file = "/ERA5_WL_SLR000_defaultdict.pkl"
wl_time_file = "/ERA5_WLtime_SLR000_defaultdict.pkl" 
tide_file = "/ERA5_tide_SLR000_defaultdict.pkl"
tide_time_file = "/ERA5_tidetime_SLR000_defaultdict.pkl"

with open(WORKING_DICTS + wl_file, 'rb') as f:
    station_data = pickle.load(f)
with open(WORKING_DICTS + wl_time_file, 'rb') as f:
    station_time = pickle.load(f)
with open(WORKING_DICTS + tide_file, 'rb') as f:
    station_data_tides = pickle.load(f)
with open(WORKING_DICTS + tide_time_file, 'rb') as f:
    station_time_tides = pickle.load(f)


############ Apply Limber Temporal Forecast of Rates given SLR #############

## for XXXcm by 2100, depending on which SLR is indicated in config file ##
cSLR = [int(SLR)]  

# get waves data
waves = xr.open_dataset(ERA5_WAVES)
wave_lons = waves['lon'].values
wave_lats = waves['lat'].values
wv_idx = []

bad_idx = np.where(np.all(np.isnan(waves['Hs'].values), axis=1))[0]

for j in range(len(df_mvar)):

    DFMlon = df_mvar['DFMlon'][j]
    DFMlat = df_mvar['DFMlat'][j]

    dist = np.sqrt((wave_lons - DFMlon)**2 + (wave_lats - DFMlat)**2)

    # exclude bad station(s)
    dist[bad_idx] = np.inf

    idx = np.argmin(dist)
    wv_idx.append(idx)
        

df_mvar['Wave_idx'] = wv_idx

## concat wave and wl data 

wl_mean_dict = {}

for key, arr_list in station_data.items():
    wl_all = np.concatenate(arr_list) / 10000
    wl_mean_dict[key] = np.nanmean(wl_all)

hs_mean = np.nanmean(waves['Hs'].values, axis=0)
tm_mean = np.nanmean(waves['Tm'].values, axis=0) #/ np.timedelta64(1, 's'), axis=0)


# define function to remove outliers
def remove_outliers_iqr(data):
    """Removes outliers from a 1D array using the IQR method."""
    q1, q3 = np.nanpercentile(data, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return data[(data >= lower_bound) & (data <= upper_bound)]

tm_err = 1.3
hs_err = 0.25
wl_err = 0.15
elev_err = 0.25
slope_err = ((elev_err * np.sqrt(2)) / df_mvar['BeachWidth_meters']) # calculate uncertainty of slope
df_mvar = df_mvar.replace('#NAME?', np.nan)
df_mvar['ActiveSlope_pct'] = pd.to_numeric( df_mvar['ActiveSlope_pct'])
beach_slope = (df_mvar['ActiveSlope_pct']/100)
df_mvar['Beach_Slope'] = beach_slope

### linear interp of MV model through time ###

def make_years(start_yr, end_yr):
    next_20 = ((start_yr + 19) // 20) * 20
    
    years = [start_yr]
    
    if next_20 != start_yr:
        years.append(next_20)
    
    years.extend(range(next_20 + 20, end_yr + 1, 20))
    
    return years

yrs = make_years(int(startyear), int(stopyear)) # each year that we will calculate a rate for, 20 yr intervals

## linear interpolation of MV regression rates 000 and xxx slr
mv_start = df_mvar[['MR','LR','UR']].values
mv_stop = df_mvar2[['MR','LR','UR']].values
mv_interp = {}
for i, yr in enumerate(yrs):
    f = f = (yr - int(startyear)) / (int(stopyear) - int(startyear))  # interpolation fraction
    mv_interp[yr] = mv_start + f * (mv_stop - mv_start)

swashes = []
setups = []

   
cSL = cSLR[0] / 100 # sea level at yr 2100, in meters
       
# slr quadratic coefficients a, b, c
A = (cSL - 0.3) / 10000
B = 0.003 - 4000 * A
C = 4000000 * A - 6

# Define Vars
g = 9.81 #m/s2, gravity
d = 1022.7 #kg/m3, density
con = (d*g)/8

df = df_mvar
df['OutputID'] = df['OutputID'].astype(int)
df = df.drop_duplicates(subset=['OutputID'], keep='first')
df.reset_index(drop=True, inplace=True)

df_w1 = pd.DataFrame() #dataframe with all vars
df_w2 = pd.DataFrame() #df with just rates and some vars

Ratevars = ['LR','MR','UR']

st = time.time()
  

j = 0
for j in range(len(df)): #for all wave file names within the table df

    DFMid = df['DFMid'][j] # wave and water level site ID
    
    # get wl stats
    wl = np.concatenate(station_data[DFMid]) / 10000
    wl = np.nanmean(wl)

    ###### DOUBLE CHECK THIS IS CORRECT CALCULATION
    MTL = np.nanmean(wl)

    idx = int(df['Wave_idx'][j])
    hs = hs_mean[idx]
    tm = tm_mean[idx]

    
    df_w1.loc[j,'Pt_ID'] = int(df.OutputID[j])
    df_w2.loc[j,'Pt_ID'] = int(df.OutputID[j])
    
    # Physical Bluff Attributes
    Toe = df.BluffToe_Elev[j]
    Beach_Slope = df.Beach_Slope[j]
    
    #TWL
    WvL = ((tm**2)*g)/(2*math.pi) # Wavelength (m), Stockdon 2006 wavelength, outputs array
    WvL_err = (tm_err * tm * g) / math.pi ## wavelength error
    setup = np.multiply((0.35*Beach_Slope),((np.multiply(hs,WvL))**0.5)) #Stockdon
    setups.append(setup)
    with np.errstate(divide='ignore', invalid='ignore'): setup_err = setup * ((slope_err[j] / Beach_Slope) + 0.5 * ((hs_err / hs) + (WvL_err / WvL)))
    setup_err = np.nan_to_num(setup_err, nan=np.nan, posinf=np.nan, neginf=np.nan) #setup error, convert invalid results (hs or WvL == 0) to NaN
    Swash = ((hs*WvL*((0.563*(Beach_Slope**2))+0.004))**0.5)/2 
    swashes.append(Swash)
    with np.errstate(divide='ignore', invalid='ignore'): swash_err = Swash * 0.5 * ((hs_err/hs) + (WvL_err/WvL) + ((1.126*slope_err[j]*Beach_Slope)/(0.563*Beach_Slope**2+0.004)))
    swash_err = np.nan_to_num(swash_err, nan=np.nan, posinf=np.nan, neginf=np.nan) #swash error, convert invalid results (hs or WvL == 0) to NaN
    
    R2 = 1.1*(setup+Swash)
    R2_err = 1.1 * np.sqrt(setup_err**2 + swash_err**2) ## R2 error
    
    ns = 0.18*Beach_Slope*tm*((g*hs)**0.5) #From Limber 2018 via Guza & Thornton 1981 for sites without beaches
    with np.errstate(divide='ignore', invalid='ignore'):ns_err = ns * ((slope_err[j]/Beach_Slope) + (tm_err/tm) + 0.5*(hs_err/hs)) # ns error
    ns_err = np.nan_to_num(ns_err, nan=np.nan, posinf=np.nan, neginf=np.nan) #ns error, convert invalid results (hs or tm == 0) to NaN
    
    # save swash and setup to dataframe
    df_w2.loc[j,'Swash_mean_meters'] = np.round(np.mean(Swash),3)
    df_w2.loc[j,'Setup_mean_meters'] = np.round(np.mean(setup),3)
    
    
    #Calibrate without SLR - Calculate Omega to calc K
    if Toe > MTL:
        #print('Toe > MTL, ' + str(df['OutputID'][j]) )
       
        TWL = (R2 + wl)/Toe #slr (m)
        TWL_err = TWL * np.sqrt((R2_err**2 + wl_err**2)/(R2 + wl)**2 + (elev_err/Toe)**2)
    else:
        print('Toe < MTL, ' + str(df['OutputID'][j]))
        TWL = (ns + wl)/Toe
        TWL_err = TWL * np.sqrt((ns_err**2 + wl_err**2)/(ns + wl)**2 + (elev_err/Toe)**2)
    
    TWL_Oi = np.nanmean(TWL)
    TWL_Oi_err = np.nanmean(TWL_err)    
    
    Hack = con*((hs+wl)**2) #(N/m)
    Hack_err = Hack * np.sqrt(hs_err**2 + wl_err**2)
    
    Hack_Oi = np.nanmean(Hack)
    Hack_Oi_err = np.nanmean(Hack_err)
    
    #Variables for Bruun Rule
    I =  df.BeachWidth_meters[j] #cross-shore length of active profile. Beach width?
    hc = df.BluffTop_Elev[j] # cliff height, or should it be cliff height from bluff toe? in m
    c = 0.5 # proportion of cliff derived sediment that remains in the nearshore 
    dep = (2.28 * hs - (68.5 * (hs**2/(g*tm**2))))  # depth of closure in m
    dep_err = np.sqrt(((2.28-(2*68.5*hs)/(g*tm**2))*hs_err)**2 + (((2*68.5*hs**2)/(g*tm**3))*tm_err)**2) # closure depth error in m
    hslr = 0.00318 # historic sea level rise rate (m/yr), or 1/8th inch per year. https://oceanservice.noaa.gov/facts/sealevel.html
    
    y = 0
    for y in range(len(Ratevars)): #For each rate (low med hi)

        rate = df.loc[j, '%s' %Ratevars[y]]
        
        #Calibrate without SLR - Calculate K
        TWL_K = rate/(TWL_Oi)
        Hack_K = rate/(Hack_Oi)

        TWL_Kerr = (rate/TWL_Oi) * (TWL_Oi_err / TWL_Oi)
        Hack_Kerr = (rate/Hack_Oi) * (Hack_Oi_err / Hack_Oi)
        
            
        df_w1.loc[j,'TWL_K'] = f'{TWL_K:.2e}'#round(TWL_K,2)
        df_w1.loc[j,'Hack_K'] = f'{Hack_K:.2e}'
        
        i = 0
        slrs = [0.0]
        slrates = [0.003]
        for i in range(len(yrs)): # For each 20-yr segment

            yr = yrs[i]
            slr = A*yr**2 + B*yr + C
            slrate = (slr - slrs[-1]) / 20 # meters / yr
            slrs.append(slr)
            slrates.append(slrate)
            
            #Total water level
            if Toe > MTL:
                TWL = (R2 + wl + (slr))/Toe #slr (m)
                TWL_err = TWL * np.sqrt((R2_err**2 + wl_err**2)/(R2 + wl)**2 + (elev_err/Toe)**2)
            else:
                TWL = (ns + wl + (slr))/Toe #From Limber 2018 via Guza & Thornton 1981 for sites without beaches
                TWL_err = TWL * np.sqrt((ns_err**2 + wl_err**2)/(ns + wl)**2 + (elev_err/Toe)**2)

            TWL_err = remove_outliers_iqr(TWL_err.flatten())
            
            #Calculate Omega
            TWL_O = np.nanmean(TWL) #Unitless
            TWL_O_err = np.nanmean(TWL_err)
            
            df_w1.loc[j,'%s_%sSLR_TWLO' %((Ratevars[y]),str(yr))] = f'{TWL_O:.2e}'#round(TWL_O,2)
    
            #Hackney
            Hack = con*((hs+wl+(slr))**2) #(N/m)
            Hack_err = Hack * np.sqrt(hs_err**2 + wl_err**2)
            
            # Calculate Omega
            Hack_O = np.nanmean(Hack)
            Hack_O_err = np.nanmean(Hack_err)
            
            df_w1.loc[j,'%s_%sSLR_HackO' %((Ratevars[y]),str(yr))] = f'{Hack_O:.2e}'
            
            # Calculate SLR and rate dependant forecasts
            TWL_KO = (TWL_K * TWL_O)
            with np.errstate(divide='ignore', invalid='ignore'): TWL_KOerr = TWL_KO * np.sqrt((TWL_Kerr/TWL_K)**2 + (TWL_O_err/TWL_O)**2)
            TWL_KOerr_err = np.nan_to_num(TWL_KOerr, nan=np.nan, posinf=np.nan, neginf=np.nan) # convert invalid results (TWL_O or TWL_K == 0) to NaN
            
            if np.isnan(TWL_KOerr):
                TWL_KOerr = 0.0
            else:           
                TWL_KOerr = TWL_KOerr 
                
            Hack_KO = (Hack_K * Hack_O)
            with np.errstate(divide='ignore', invalid='ignore'): Hack_KOerr = Hack_KO * np.sqrt((Hack_Kerr/Hack_K)**2 + (Hack_O_err/Hack_O)**2)
            Hack_KOerr = np.nan_to_num(Hack_KOerr, nan=np.nan, posinf=np.nan, neginf=np.nan) # convert invalid results (Hack_K or Hack_O == 0) to NaN
            
            if np.isnan(Hack_KOerr):
                Hack_KOerr = 0.0
            else:            
                Hack_KOerr = Hack_KOerr

            # calculate bruun rates 
            dSLR =  slrate - 0.003
            bruun_rates = rate + ((dSLR * I) / (dep + c*hc)) * 100 # multiply by 100 for m to cm conversion
            I_err = (np.sqrt(elev_err**2 + dep_err**2)) / (df.Beach_Slope[j]) 
            bruun_errs = np.sqrt(((dSLR*I_err)/(c*hc+dep))**2 + ((dSLR*I*c*elev_err)/(c*hc+dep)**2)**2 + ((dSLR*I*dep_err)/(c*hc+dep)**2)**2)
            bruun_err = np.nanmean(bruun_errs)
            bruun_rate = np.nanmean(bruun_rates)
            if bruun_rate < 0.0:
                bruun_rate ==0
            else:
                bruun_rate = bruun_rate
            
            ## get MV model values for weighting
            mv_vals = mv_interp[yr][j]
            MV_MR = mv_vals[0]
            MV_LR = mv_vals[1]
            MV_UR = mv_vals[2]
            if Ratevars[y] == 'LR':
                MV_rate = MV_LR
            elif Ratevars[y] == 'MR':
                MV_rate = MV_MR
            elif Ratevars[y] == 'UR':
                MV_rate = MV_UR
            MV_err = (MV_UR - MV_LR) / 2    
                
            ## Calculate weighting for each model
            IcM = np.abs(MV_rate - ((Hack_KO + TWL_KO + bruun_rate) / 3))
            IcT = np.abs(TWL_KO - ((Hack_KO + bruun_rate + MV_rate)/3))
            IcH = np.abs(Hack_KO - ((TWL_KO + bruun_rate + MV_rate)/3))
            IcB = np.abs(bruun_rate - ((Hack_KO + TWL_KO + MV_rate)/3))
            yM = 1 / (IcM + 2*MV_err)**2
            yT = 1 / (IcT + 2*TWL_KOerr)**2
            yH = 1 / (IcH + 2*Hack_KOerr)**2
            yB = 1 / (IcB + 2*bruun_err)**2

            Final_Pred = (
                Hack_KO*yH +
                TWL_KO*yT +
                bruun_rate*yB +
                MV_rate*yM
            ) / (yH + yT + yB + yM)
            
            Final_Uncert = 1 / np.sqrt(yH + yT + yB + yM)
            
            ## write final predictions to dataframe

            if Ratevars[y] == 'LR':

                df_w1.loc[j,'%s_%sSLR_TWL' %((Ratevars[y]),(str(yr)))] = np.around((TWL_KO - TWL_KOerr),2)
                df_w1.loc[j,'%s_%sSLR_Hack' %((Ratevars[y]),(str(yr)))] = np.around((Hack_KO - Hack_KOerr),2)
                df_w1.loc[j,'%s_%sSLR_Bruun' %((Ratevars[y]),(str(yr)))] = np.around((bruun_rate - bruun_err),2)
                df_w1.loc[j,'%s_%sSLR_MV' %((Ratevars[y]),(str(yr)))] = np.around((MV_rate - MV_err),2)

                if Final_Uncert > Final_Pred:
                    Final_Uncert = Final_Pred
                else:
                    Final_Uncert = Final_Uncert
                df_w1.loc[j,'%s_%sSLR_Rate' %((Ratevars[y]),(str(yr)))] = np.around((Final_Pred - Final_Uncert),2)
                df_w2.loc[j,'%s_%sSLR_Rate' %((Ratevars[y]),(str(yr)))] = np.around((Final_Pred - Final_Uncert),2)
            
            elif Ratevars[y] == 'MR':

                df_w1.loc[j,'%s_%sSLR_TWL' %((Ratevars[y]),(str(yr)))] = np.around((TWL_KO),2)
                df_w1.loc[j,'%s_%sSLR_Hack' %((Ratevars[y]),(str(yr)))] = np.around((Hack_KO),2)
                df_w1.loc[j,'%s_%sSLR_Bruun' %((Ratevars[y]),(str(yr)))] = np.around((bruun_rate),2)
                df_w1.loc[j,'%s_%sSLR_MV' %((Ratevars[y]),(str(yr)))] = np.around((MV_rate),2)
                
                df_w1.loc[j,'%s_%sSLR_Rate' %((Ratevars[y]),(str(yr)))] = np.around((Final_Pred),2)
                df_w2.loc[j,'%s_%sSLR_Rate' %((Ratevars[y]),(str(yr)))] = np.around((Final_Pred),2)
                
            
            elif Ratevars[y] == 'UR':

                df_w1.loc[j,'%s_%sSLR_TWL' %((Ratevars[y]),(str(yr)))] = np.around((TWL_KO + TWL_KOerr),2)
                df_w1.loc[j,'%s_%sSLR_Hack' %((Ratevars[y]),(str(yr)))] = np.around((Hack_KO + Hack_KOerr),2)
                df_w1.loc[j,'%s_%sSLR_Bruun' %((Ratevars[y]),(str(yr)))] = np.around((bruun_rate + bruun_err),2)
                df_w1.loc[j,'%s_%sSLR_MV' %((Ratevars[y]),(str(yr)))] = np.around((MV_rate + MV_err),2)

                df_w1.loc[j,'%s_%sSLR_Rate' %((Ratevars[y]),(str(yr)))] = np.around((Final_Pred + Final_Uncert),2)
                df_w2.loc[j,'%s_%sSLR_Rate' %((Ratevars[y]),(str(yr)))] = np.around((Final_Pred + Final_Uncert),2)
    
    et = time.time()
    sys.stdout.write('\rElapsed time: %.ds; Sites Completed: %.d; Remaining Sites: %s' %((et-st),(j+1),(len(df)-(j+1))))
    sys.stdout.flush()

df_w1['OutputID'] = df_w1['Pt_ID'].astype(int)
df_w2['OutputID'] = df_w2['Pt_ID'].astype(int)


### Save Dataframe ####
fdf = os.path.join(WORKING_DIR,'{}.csv'.format(dfN_Fcast)) #saves dataframe to specified folder
fdf2 = os.path.join(WORKING_DIR,'{}.csv'.format(dfN_Fcast2)) #saves dataframe to specified folder
df_Fcast = pd.merge(df, df_w1, on="OutputID", how='left') #merge with intial dataframe, match by map label
df_Fcast2 = pd.merge(df, df_w2, on="OutputID", how='left') #merge with intial dataframe, match by map label
df_Fcast.to_csv(str(fdf))
df_Fcast2.to_csv(str(fdf2))

##### Save the rates dataframe as a shapefile ####

geometry = [Point(xy) for xy in zip(df_Fcast2['BluffTop_UTMX'], df_Fcast2['BluffTop_UTMY'])]
# Create a GeoDataFrame
gdf = gpd.GeoDataFrame(df_Fcast2, geometry=geometry)
# Set the CRS to UTM Zone 10N, NAD83(2011), NAVD88
gdf.crs = "EPSG:6339"
# Save to shapefile
fShp2 = os.path.join(WORKING_DIR,'{}.shp'.format(dfN_Fcast2)) #saves shapefile to specified folder
gdf.to_file(str(fShp2))
