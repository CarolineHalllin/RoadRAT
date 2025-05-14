# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 13:57:09 2023

@author: Caroline Hallin

Script to define dictionaries to store all model parameters.
"""

# Create empty dictionaries

# Keys in dictionary for road point information:
s = {
    'ID': [],           # Index road point [integer]
    'x': [],            # X-coordinate road point, m [Float]
    'y': [],            # Y-coordinate road point, m [Float]
    'z': [],            # Z-coordinate road point, m rel elevation reference [Float]
    'Lveg': [],         # Distance from road to vegetation line, m [Float]
    'Lveg_hist': [],    # Distance from road to historical vegetation line, m [Float]
    'Lveg_slr': [],     # Distance from road to vegetation line in future scenario, m [Float]
    'B': [],            # Beach width defined as distance from veg line to shoreline m [Float]
    'V': [],            # Volume of sediment between the road and vegetation line, m3 [Float]
    'V_r': [],          # Volume of sediment between the road and 10-year runup level, m3 [Float]
    'V_slr': [],        # Volume of sediment between the road and vegetation line in future scenario, m3 [Float]
    'V_r_slr': [],      # Volume of sediment between the road and 10-year runup limit line in future scenario, m3 [Float]
    'h': [],            # Average height between vegetation line and road, rel vegetation line [Float]
    'h_slr': [],        # Average height between vegetation line and road in future scenario, rel vegetation line [Float]
    'h_max': [],        # Maximum height between vegetation line and road, rel elevation of vegetation line [Float]
    'h_max_slr': [],    # Maximum height between vegetation line and road in future scenario, rel elevation of vegetation line [Float]
    'alpha': [],        # Slope between shoreline and the vegetation line, m/m [Float]
    'sh_dir': [],       # Orientation of shore normal, degrees [Float]
    'x_wave':[],        # X-coordinate wave data extraction point, m [Float]
    'y_wave':[],        # Y-coordinate wave data extraction point, m [Float]
    'D50': [],          # Median grain size in material eroded by, m (based on soil class)
    'c_impact': [],     # Erosion coefficient (based on soil class ) Larson impact formula, -, [Float] 
    'struc': [],        # Non erodible layer between road and shoreline y/n [Bool]
    'Dc': [],           # Closure depth, m [Float]
    'Db': [],           # Dune toe elevation (elevation of vegetation line), m rel ref [Float]
    'Db_r': [],         # Dune toe elevation based on average yearly runup level, m rel ref [Float]
    'LDean': [],        # Length of Dean profile, m [Float]
    'p_swl': [],        # Probability of flooding due to SWL today, return period years [Float]
    'p_r': [],          # Probability of flooding due to runup today, return period years [Float]
    'p_er': [],         # Probability of erosion today, return period years [Float]
    'p_swl_slr': [],    # Probability of flooding due to SWL in the future, return period years [Float]
    'p_r_slr': [],      # Probability of flooding due to runup in the future, return period years [Float]
    'p_er_slr': [],     # Probability of erosion in the future, return period years [Float]
    'swl_100': [],      # 100-year SWL, m [Float]
    'r_100': [],        # 100-year TWL, m [Float]
    'er_100': [],       # 100-year erosion, m3/m [Float]
    'er_slr_100': []   # 100-year erosion in the future, m3/m [Float]
}

# Keys in dictionary for values that are constant for all road points and either defined by user or calculated in the program:
c = {
    'spacing': None,    # Road point spacing, m [Float]
    'start_yr': None,   # Start year for analysis [integer]
    'end_yr': None,     # End year for analysis [integer]
    'SLR': None,        # Relative SLR until endyear, m [Float]
    'MSL': None,        # Mean sea level, m relative elevation ref [Float]
    'FS_depth': None,   # Depth of deepest part of the foreshore (where FS_file is extracted), m rel MSL [Float]
    'Db_max' : None,    # Maximum dune toe elevation, m rel elevation reference [Float]
    'CL_yr': None,      # Year of historical coastline [integer]
    'depth_wave_file': None, # Depth at wave data extraction point, m rel elevation reference [Float]
    'output_file': None, # Output file name
    'crs' : 'EPSG:3006', # crs code for the pandas geodataframe, default SWEREF 99 TM
    'crs_wave' : 'EPSG:4326', #crs code for wave input data, default WGS84
    'd50_default': 0.2*10**-3, # Default median grain size, m [Float]
    'c_impact_default': 0.0003,   # Default erosion coefficient, - [Float]
    'viscosity': 1.3*10**-6, # Kinematic viscosity of water, m2/s [Float]
    'SWL_1_x': None,    # X-coordinate SWL point 1, m [Float]
    'SWL_1_y': None,    # Y-coordinate SWL point 1, m [Float]
    'SWL_2_x': None,    # X-coordinate SWL point 2, m [Float]
    'SWL_2_y': None,    # Y-coordinate SWL point 2, m [Float]
    'road_width': None # Width of road, m [Float]
}

# Keys in dictionary with input file names:
f = {
    'SWL_file_1': None,   # Textfile with observed water levels [time, WL(m)]
    'SWL_file_2': None,   # Textfile with observed water levels [time, WL(m)]
    'road_file':None,  # Shp-file with road line segments (most nearshore road location)
    'DEM_file': None,   # Digital elevation model extending from the shoreline to the road
    'CL_file': None,    # Shapefile with a line indicating the vegetation line
    'SH_file': None,    # Shapefile with a line indicating the shoreline
    'FS_file': None,    # Shapefile with a line indicating the end of the foreshore
    'hist_CL_file': None, # Shapefile with a line indicating the historical vegetation line
    'NE_file': None, # Shapefile with line segments indicating revetments
    'grain_size_file': None,  # Shapefile with polygons - d50 specification for Dean profile calculations
    'erodibility_file': None, # Shapefile with polygons - dpecification of Cimpact coefficient for storm erosion
    'wave_files_dir': None,  # Directory with Netcdf file with wave data
    'wave_file': None, #File with wave data
    'bathyline_file': None # Shapefile with bathyline indicating depth contour for wave data extraction
}

# Keys in dictionary for process information:
p = {
    'process_inundation': True,  # Process inundation, True or False [Bool]
    'process_erosion': True,    # Process erosion, True or False [Bool] (if True, process_runup is also True)
    'process_runup': True,      # Process runup, True or False [Bool]
    'process_SLR': True,         # Process SLR, True or False [Bool]
    'process_extract_waves': True, #Indicate if waves should be extracted from netcdf-files within RoadRAT: else provided in file
    'process_plotting': True,    # Process plotting, True or False [Bool]
    'process_calibration': False, # Process calibration, True or False [Bool]
    'process_hallermeier_rhodes': True # Process Hallermeier and Rhodes, True or False [Bool]
}
