### final code, compute the setback distance of future bluff crest based on changing rate through time
### manually adjust the year intervals (lines 52 - 60) if you used a different start year and end year than the example defaults

import geopandas as gpd
from shapely.geometry import Point
import math
from collections import defaultdict
import numpy as np

from config import SLR, WORKING_DIR, FINAL_DIR, county, stopyear

# Function to calculate new point based on bearing and distance in UTM coordinates
def calculate_new_point(x, y, bearing, distance):
    # Convert bearing to radians
    bearing = math.radians(bearing)
    # Calculate the new coordinates
    new_x = x + distance * math.sin(bearing)
    new_y = y + distance * math.cos(bearing)
    return Point(new_x, new_y)
    
# Load the shapefile
path = WORKING_DIR + '/JustEnsembleRates_' + county + '_' + SLR + 'cm_by' + stopyear +'.shp'
crest = gpd.read_file(path)

# Define rate labels and append
ratelabel = ['LR','MR', 'UR']
append = ['_Lower','', '_Upper']

distancedict = defaultdict(list)


# Main loop for processing
for i in range(0, 3):  # Loop over three labels
    print(str(i), flush=True)
    # Create a list for the new points
    new_points = []
    dist = []
    
    # Iterate through rows in the GeoDataFrame
    for index, row in crest.iterrows():
        # Extract the original geometry and attributes
        original_point = row.geometry
        
        bearing = row['Bearing']
        label = ratelabel[i]
        
        # Get rates from dataframe
        t017 = row[label + '_2017SLR'] / 100
        t020 = row[label + '_2020SLR'] / 100
        t040 = row[label + '_2040SLR'] / 100
        t060 = row[label + '_2060SLR'] / 100
        t080 = row[label + '_2080SLR'] / 100
        t100 = row[label + '_2100SLR'] / 100
        
        # calculate the setback distance 
        setback = t017 * 3 + t020 * 13 + t040 * 20 + t060 * 20 + t080 * 20 + t100 * 10 

        dist.append(setback)
    
        # Calculate the new point
        new_point = calculate_new_point(original_point.x, original_point.y, bearing, -setback)
        new_points.append(new_point)
        
    distancedict[ratelabel[i]].append(dist)
    
    # Create a new GeoDataFrame with updated points
    new_gdf = gpd.GeoDataFrame(geometry=new_points, crs=crest.crs)
    
    # Add only the columns you want to keep
    #new_gdf['TransectID'] = crest['TransectID']   # <-- keep this column
    new_gdf['OutputID'] = crest['TransectID']     # <-- keep this column
    new_gdf['SectionID'] = crest['SectionID']         # <-- keep this column
    new_gdf['distance'] =  np.round(dist,2)
    
    # Save the resulting shapefile. Separate shapefiles for best estimate, upper uncertainty range, and lower uncertainty range
    output_path = FINAL_DIR + '/' + county + '_SLR' + SLR + 'cm_by' + stopyear + append[i] + '_4ModelEnsemble.shp'
    new_gdf.to_file(output_path)