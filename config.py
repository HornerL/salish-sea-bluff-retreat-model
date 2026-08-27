
from pathlib import Path

### Project settings

SLR = "100" # projected sea level rise by 2100, in cm

### Input data

# file containing known erosion rates
FIELD_DATA = Path(r"C:\Users\lhorner\Data\Bluff_erosion_CGS\Snohomish_forecasted_20260420\Data\ValidationSites_BoundaryConditions_20260417.csv")

# transects for erosion rate predictions
TRANSECT_DATA = Path(r"C:\Users\lhorner\Data\Bluff_erosion_CGS\IslandCo_forecasted_20260612\Data\IslandCo_BoundConditions_20260612.csv")

# CoSMoS wave conditions, 000 m SLR
ERA5_WAVES = Path(r"F:\CoSMoS2026\Waves\Island\Reanalysis_and_Projected_CoSMoSwaves_IslandCounty_sealevel000m.nc") 

# CoSMoS waterelevel conditions, 000 m SLR
ERA5_WL_DIR = Path(r"F:\CoSMoS2026\WaterLevels\ERA5\ERA5_000")

# CoSMoS tides conditions, 000 m SLR
ERA5_TIDE_DIR = Path(r"F:\CoSMoS2026\Tides\ERA5_tidal_Results\000")

# Digital Elevation Model (USGS CoNED)
CONED_DTM = Path(r"E:\USGS_laptop\CoNED\new_CoNED_king_pierce\tSound_CoNED_Edit9.tif")

# directory where computed boundary conditions will be saved
WORKING_DIR = Path(r"C:\Users\lhorner\Data\Bluff_erosion_CGS\IslandCo_forecasted_20260612\working_data")

# filename that the computed boundary conditions will be saved to
dfN_i = 'IslandCo_Computed_BoundConds_2026xxxx'