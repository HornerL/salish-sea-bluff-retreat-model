
### Project settings

# county or study area name
county = "Skagit"

# year that initial bluff conditions are defined (based on the LiDAR)
startyear = '2017'

# year for future SLR projections
stopyear = '2100'

SLR = "100" # projected sea level rise in cm, by stopyear

### Input data and directories

# file containing known erosion rates
FIELD_DATA = r"C:\Users\lhorner\Data\Bluff_erosion_CGS\Snohomish_forecasted_20260420\Data\ValidationSites_BoundaryConditions_20260417.csv"

# transects for erosion rate predictions
TRANSECT_DATA = r"C:\Users\lhorner\Data\Bluff_erosion_CGS\Skagit_forecasted_20260902\Data\Skagit_BoundConds_20260902.csv"

# CoSMoS wave conditions, 000 m SLR
ERA5_WAVES = r"E:\CoSMoS2026\Waves\Skagit\Reanalysis_and_Projected_CoSMoSwaves_Skagit_sealevel000m.nc"

# CoSMoS waterelevel conditions, 000 m SLR
ERA5_WL_DIR = r"E:\CoSMoS2026\WaterLevels\ERA5\ERA5_000"

# CoSMoS tides conditions, 000 m SLR
ERA5_TIDE_DIR = r"E:\CoSMoS2026\Tides\ERA5_tidal_Results\000"

# Digital Elevation Model (USGS CoNED)
CONED_DTM = r"D:\USGS_laptop\CoNED\new_CoNED_king_pierce\tSound_CoNED_Edit9.tif"

# directory where computed intermediate dataframes will be save
WORKING_DIR = r"C:\Users\lhorner\Data\Bluff_erosion_CGS\Skagit_forecasted_20260902\working_data"

# directory where final dataframes will be saved
FINAL_DIR = r"C:\Users\lhorner\Data\Bluff_erosion_CGS\Skagit_forecasted_20260902\Data"

# filename that the computed boundary conditions will be saved to
dfN_i = 'SkagitCo_ERA5_Computed_BoundConds_20260902'

# filename for 1st round of historical rate predictions (000 m SLR)
hist_rates_file = 'SkagitCo_ERA5_HistoricalRates_20260902'

# directory where compiled time series wl/tide data will be stored (large files)
WORKING_DICTS = r"D:\USGS_laptop\CoSMoS_WaveWL_outputs\working_dictionaries\Skagit"

# directory of CMIP6 future waves file
CMIP6WAVES = r"E:\CoSMoS2026\Waves\Skagit\Reanalysis_and_Projected_CoSMoSwaves_Skagit_sealevel"

# directory of CMIP6 future wl file
CMIP6WL_DIR = r"E:\CoSMoS2026\WaterLevels\ERA5\ERA5_"

# directory of CMIP6 future tides file
CMIP6TIDE_DIR = r'E:\CoSMoS2026\Tides\ERA5_tidal_Results'

# directory of cmip_diff file
CMIPDIFF_DIR = r"E:\CoSMoS2026\WaterLevels\CMIP6_cdf_diff"