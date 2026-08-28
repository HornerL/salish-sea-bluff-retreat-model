import pandas as pd
import numpy as np
import time
import sys
import json
import statsmodels.formula.api as smf

from config import dfN_i, WORKING_DIR, FIELD_DATA, hist_rates_file

df = pd.read_csv(WORKING_DIR + '/' + dfN_i + '.csv')

df_field = pd.read_csv(FIELD_DATA)

df_field['NTR95_per_TWLrange'] = df_field['NTR95'] / df_field['TWL_Range']


### Compute multivariate OLS regression from field data, initially without bootstrapping ###

model = smf.ols(formula = 'ErosionRate_observed ~ NTR95_per_TWLrange + TWL_Skew + WavePower_Cum +  WavePower_Skew + Cohesion_MidPoint + Cohesion_Range', data=df_field)
res = model.fit()
print(res.summary())

## drop NaN rows from df_Vars (ERIC TO CHECK COUNTY TRANSECT FOR NANS)
df_Vars = df
predictor_columns = [ 'NTR95_per_TWLrange', 'TWL_Skew', 'WavePower_Cum',  'WavePower_Skew', 'Cohesion_MidPoint', 'Cohesion_Range']
df_Vars = df_Vars.dropna(subset=predictor_columns)

## create dataframe of the non-bootstrapped rate predictions, for plotting and comparison
pred = res.get_prediction(df_Vars)
nonstrap = pred.summary_frame()
nonstrap['TransectID'] = df_Vars['TransectID']


# ## Now we apply the regression to all our transects
# We perturb each variable by its known measurment error and run the regression for 1000 iterations, then calculate the mean EPR for each transect
# For upper and lower uncertainty, we find the 97.5th and 2.5th percentile of the prediction intervals from the 1000 iterations

## bootsrapped rate predictions, and prediction intervals ##

## additional measurement error/uncertainty
EPRerror_DSAS = 0 # assume negligable uncertainty for field-measured rates. OR 2.0 cm/yr for DSAS-derived rates
wl_error = 0.3 # meters (wl error is 0.15 m, but we are propogating it across range from min to max)
co_error = 23 # cohesive strength error
NTRperWL_error = np.nanmean(df_field['NTR95_per_TWLrange']) * np.sqrt((wl_error/df_field['NTR95'])**2 + 
                                                                      (wl_error/df_field['TWL_Range'])**2)
with open(WORKING_DIR + '/wp_uncertainty.json', "r") as f:
    parameters = json.load(f)
wp_error = parameters["relative_wave_power_uncertainty"]

## vary values by uncertainties 

def vary(series, uncertainty):
    ##Varies values in a pandas Series by a random amount up to max_variation.
    return series + np.random.normal(loc=0, scale=uncertainty/3, size=len(series))

pred_res = {}
res_list = []
r2_list = []
r2adj_list = []
pval_list = []
aic_list = []
pred_list = []
piu_list = []
pil_list = []

## bootstraping to predict new transects with uncertainty

df_sort = df_field.sort_values(by='ErosionRate_observed', ascending=True).reset_index(drop=True) 
st = time.time()

for i in range (0, 1000):
    st1 = time.time()
    df_pert = df_sort.copy()
    
    ## vary these values by +/- uncertainty
    EPR = vary(df_pert['ErosionRate_observed'], EPRerror_DSAS) 
    WvPwr = vary(df_pert['WavePower_Cum'], (df_pert['WavePower_Cum'] * wp_error))
    NTRperWL = vary(df_pert['NTR95_per_TWLrange'], NTRperWL_error)
    cohes = vary(df_pert['Cohesion_MidPoint'], co_error)
        
    ## make sure none of the wave power and EPR values are less than 0
    NTRperWL[NTRperWL < 0] = 0.0
    df_pert['NTR95_per_TWLrange'] = NTRperWL

    WvPwr[WvPwr < 0] = 0.0
    df_pert['Cumul_WavePower'] = WvPwr

    EPR[EPR < 0] = 0.0
    df_pert['ErosionRate_observed'] = EPR

    cohes[cohes < 0] = 0.0
    df_pert['Cohesion_Mean'] = cohes

    ## fit regression to new adjusted data    
    model = smf.ols(formula = 'ErosionRate_observed ~ NTR95_per_TWLrange + TWL_Skew + WavePower_Cum +  WavePower_Skew + Cohesion_MidPoint + Cohesion_Range' , data=df_pert)
    res = model.fit()

    r2_list.append(res.rsquared)
    r2adj_list.append(res.rsquared_adj)
    pval_list.append(res.f_pvalue)
    aic_list.append(res.aic)
    
    pred = res.get_prediction(df_Vars) # predict the new transect rates
    #pred = res.get_prediction(df_pert)
    pred_res[i] = pred.summary_frame()
    res_list.append(pred_res[i])
    pil_list.append(pred_res[i]['obs_ci_lower'].to_list())
    piu_list.append(pred_res[i]['obs_ci_upper'].to_list())
    pred_list.append(pred_res[i]['mean'].to_list())

    et = time.time()
    sys.stdout.write('\rElapsed time: %.ds; Estimated time remaining: %.2ds; Iterations Completed: %.d; Remaining Iterations: %s' %((et-st),((et-st1)*(1000-(i+1))),(i+1),(1000-(i+1))))
    sys.stdout.flush()
    

## R-squared of bootstrapped regression
print('R-squared:' + str(np.round(np.mean(r2_list), 3)))


## get 95th percent interval of bootstrapped prediction intervals

pil_25 = []
piu_975 = []

for i in range (0, len(pil_list[0])):
    indlist_l = []
    indlist_u = []
    for j in range (0,1000):
        low = pil_list[j][i]
        high = piu_list[j][i]
        indlist_l.append(low)
        indlist_u.append(high)
    pctl_l = np.percentile(indlist_l, 2.5)
    pctl_u = np.percentile(indlist_u, 97.5)
    pil_25.append(pctl_l)
    piu_975.append(pctl_u)
    
## Concatenate all DataFrames along a new axis
combined_df = pd.concat(res_list)

## Group by index and compute the mean for each cell
averaged_res = combined_df.groupby(combined_df.index).mean()

## combine bootstrapped values to single dataframe

bs_processed = averaged_res['mean'].to_frame(name = 'MR')
bs_processed['LR'] = pil_25
bs_processed['UR'] = piu_975
bs_processed['TransectID'] = df_Vars['TransectID']

## merge bootstrapped predictions to main variables dataframe
df_mvar = pd.merge(bs_processed, df, on='TransectID', how='right')

## convert negative rates to 0
df_mvar['LR'] = np.where(df_mvar['LR'] < 0, 0.0, df_mvar['LR'])

## save historic bootstrap projections
df_mvar.to_csv(WORKING_DIR + '/' + hist_rates_file + '.csv', index=False)