
### Project settings

SLR = "100" # projected sea level rise by 2100, in cm

### Input data

# file containing known erosion rates
FIELD_DATA = r"C:\Users\lhorner\Data\Bluff_erosion_CGS\Snohomish_forecasted_20260420\Data\ValidationSites_BoundaryConditions_20260417.csv"

# transects for erosion rate predictions
TRANSECT_DATA = r"C:\Users\lhorner\Data\Bluff_erosion_CGS\IslandCo_forecasted_20260612\Data\IslandCo_BoundConditions_20260612.csv"

# CoSMoS wave conditions, 000 m SLR
ERA5_WAVES = r"F:\CoSMoS2026\Waves\Island\Reanalysis_and_Projected_CoSMoSwaves_IslandCounty_sealevel000m.nc"

# CoSMoS waterelevel conditions, 000 m SLR
ERA5_WL_DIR = r"F:\CoSMoS2026\WaterLevels\ERA5\ERA5_000"

# CoSMoS tides conditions, 000 m SLR
ERA5_TIDE_DIR = r"F:\CoSMoS2026\Tides\ERA5_tidal_Results\000"

# Digital Elevation Model (USGS CoNED)
CONED_DTM = r"E:\USGS_laptop\CoNED\new_CoNED_king_pierce\tSound_CoNED_Edit9.tif"

# directory where computed intermediate dataframes will be save
WORKING_DIR = r"C:\Users\lhorner\Data\Bluff_erosion_CGS\IslandCo_forecasted_20260612\working_data"

# directory where final dataframes will be saved
FINAL_DIR = r"C:\Users\lhorner\Data\Bluff_erosion_CGS\IslandCo_forecasted_20260612\Data"

# filename that the computed boundary conditions will be saved to
dfN_i = 'IslandCo_ERA5_Computed_BoundConds_20260612'

# filename for 1st round of historical rate predictions (000 m SLR)
hist_rates_file = 'IslandCo_ERA5_HistoricalRates_20260612'

# directory where compiled time series wl/tide data will be stored (large files)
WORKING_DICTS = r"E:\USGS_laptop\CoSMoS_WaveWL_outputs\working_dictionaries\Island"

# directory of CMIP6 future waves file
CMIP6WAVES = r"F:\CoSMoS2026\Waves\Island\Reanalysis_and_Projected_CoSMoSwaves_IslandCounty_sealevel"

# directory of CMIP6 future wl file
CMIP6WL_DIR = r"F:\CoSMoS2026\WaterLevels\ERA5\ERA5_"

# directory of CMIP6 future tides file
CMIP6TIDE_DIR = r'F:\CoSMoS2026\Tides\ERA5_tidal_Results'

# directory of cmip_diff file
CMIPDIFF_DIR = r"F:\CoSMoS2026\WaterLevels\CMIP6_cdf_diff"