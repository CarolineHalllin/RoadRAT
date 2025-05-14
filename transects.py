# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 13:57:09 2023

@author: Caroline Hallin

Script to define functions to extract computation points along the road and calculate morphological parameters.
"""

import numpy as np
import pandas as pd
import shapefile
import rasterio
import matplotlib.pyplot as plt
import geopandas as gpd
import time
from rasterio.windows import Window
from rasterio.plot import show
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points
from input import *
from dictionaries import *
from longterm import *


# Generate computation points along the road and extract x, y and z coordinates
def computation_points(s,c,f):
# Open the shapefile
    sf = shapefile.Reader(f['road_file'])

    # Initialize lists to store results for all lines
    all_lines = []
    all_total_lengths = []
    all_distances = []
    all_ids = []
    all_x_coords = []
    all_y_coords = []

    # Initialize a global counter for unique IDs
    global_id_counter = 1

    # Iterate over all shapes in the shapefile
    for shape_index, shape in enumerate(sf.shapes()):
        points = shape.points
        
        # Create a LineString object from the points
        line = LineString(points)
        all_lines.append(line)

        # Calculate the total length of the line
        total_length = line.length
        all_total_lengths.append(total_length)

        # Generate points every x meter as defined by c['spacing']
        spacing = int(c['spacing'])
        distances = [i for i in range(0, int(total_length), spacing)]
        # Handle the case when the line is shorter than the spacing
        if not distances:
            distances = [0, total_length]  # Add start and end points of the line



        all_distances.append(distances)

        # Create unique IDs for the computation points
        ids = [global_id_counter + i for i in range(len(distances))]
        global_id_counter += len(distances)
        all_ids.extend(ids)

        # Extract x and y coordinates at regular intervals
        extracted_points = [line.interpolate(distance) for distance in distances]
        x_coords, y_coords = zip(*[(point.x, point.y) for point in extracted_points])
        all_x_coords.extend(x_coords)
        all_y_coords.extend(y_coords)

    # Store the results in the dictionary `s`
    s['ID'] = all_ids
    s['x'] = all_x_coords
    s['y'] = all_y_coords

    # Extract z coordinates from the DEM
    # Read the DEM file
    with rasterio.open(f['DEM_file']) as dem:
        # Convert x, y coordinates to row, col indices and extract z-values
        rowcol_pairs = [dem.index(x, y) for x, y in zip(s['x'], s['y'])]
        rows, cols = zip(*rowcol_pairs)

        dem_data = dem.read(1)
        z_values = [float(dem_data[row, col]) for row, col in zip(rows, cols)]
        s['z'] = z_values

    if p['process_plotting']:
        # Plot the road and computation points
        plot_computation_points(s,c,f)

        #Calculate the coordinates on the bathyline for wave extraction
    if (p['process_SLR'] or p['process_runup'] or p['process_erosion']) and p['process_extract_waves']:
        bathyline_coords(s,c,f)

    #Check if there is a non-erodible feature between the road and the shoreline
    if f['NE_file'] is not None and (p['process_runup'] or p['process_erosion']):
        check_struc(s,c,f)
    else:
        s['struc'] = [False]*len(s['x'])

        
    return s

# Extract morphological parameters from shapefiles    
def morph_parameters(s,c,f,p):
    #Generate transects, vegetation line and shoreline
    transects, cl_line, sh_line, fs_line = generate_transects(f)
    
    # Initialize the lists to store morphological parameters
    s['Lveg'] = [None]*len(s['x'])
    s['B'] = [None]*len(s['x'])
    s['V'] = [None]*len(s['x'])
    s['h'] = [None]*len(s['x'])
    s['h_max'] = [None]*len(s['x'])
    s['alpha'] = [None]*len(s['x'])
    s['Db'] = [None]*len(s['x'])
    s['sh_dir'] = [None]*len(s['x'])


    # Iterate through the road computation points
    for index in range(len(s['x'])):
        if s['struc'][index]:
            s['Lveg'][index] = None
            s['B'][index] = None
            s['V'][index] = None
            s['h'][index] = None
            s['h_max'][index] = None
            s['alpha'][index] = None
            s['Db'][index] = None
        else:
        
            road_point = Point(s['x'][index], s['y'][index])
            
            # Find the corresponding transect for the road point
            transect = transects.iloc[index].geometry
            elevations = transects.iloc[index].elevations
            
            # Calculate the intersection point between the transect and the vegetation line
            intersection = transect.intersection(cl_line)
            #Coordinates of the intersection between the transect and the shoreline
            intersection_sh = Point(transect.coords[-1])

            # Calculate the orientation of a line normal to the shoreline in the intersection point
            
            # Snip the intersection point to the shoreline

            nearest_shore_point = nearest_points(intersection_sh, sh_line)[1]
            
            # Generate two points on the shoreline at 10 m distance to each side of the nearest shore point
            distances = [-10, 10]
            points = [sh_line.interpolate(sh_line.project(nearest_shore_point) + d) for d in distances]
            
            # Draw a line between these two points
            line = LineString(points)
            
            # Rotate the line 90 degrees to get a line normal to the shoreline
            midpoint = line.interpolate(0.5, normalized=True)
            dx = points[1].x - points[0].x
            dy = points[1].y - points[0].y
            angle = np.degrees(np.arctan2(dy, dx)) + 90
            angle_rad = np.radians(angle)

            rotated_point = Point(
                midpoint.x + 10 * np.cos(angle_rad),
                midpoint.y + 10 * np.sin(angle_rad)
            )

            rotated_line = LineString([midpoint, rotated_point])

            # Measure the orientation of the rotated line in degrees north
            orientation = calculate_orientation(rotated_line)

            # Modify the orientation to be in the direction away from the road point
            if rotated_line.distance(road_point) < line.distance(road_point):
                orientation = (orientation + 180) % 360 # Add 180 degrees and take the modulo 360 to get the orientation away from the road point

            s['sh_dir'][index] = orientation


            # Calculate the distance from the road point to the intersection
            if not intersection.is_empty and intersection.geom_type == 'Point':
                distance = road_point.distance(intersection)
            else:
                distance = None  # Handle cases where there is no intersection
            
            #Calculate the distance from the vegetation line to the shoreline
            beach_width = intersection.distance(intersection_sh)

            # Calculate the volume of sediment between the road and the vegetation line
            if distance is not None:     
                # Calculate the intersection elevation
                intersection_elevation = elevations[int(distance * 2)]

                # Limit the intersection elevation to Db_max
                intersection_elevation = min(intersection_elevation, float(c['Db_max']))
        
                # Filter the elevations between the road and vegetation lines and subtract the intersection elevation from the elevations
                start_index = int(c['road_width']) # Start half the road distance from the road point (Divided by two to get half the road width and multiplied by two because grid size i 0.5 m)
                filtered_elevations = [elev - intersection_elevation for elev in elevations[start_index:int(distance * 2)] if elev - intersection_elevation >= 0]

                if filtered_elevations:
                    # Calculate the average height of the filtered elevations
                    h = np.mean(filtered_elevations)
                    # Calculate the maximum height of the filtered elevations
                    h_max = max(filtered_elevations)
                    # Calculate the volume of sediment above the intersection elevation
                    V = np.trapz(filtered_elevations, dx=0.5)
                else:
                    h = h_max = V = 0
        
            # Calculate the slope between the shoreline and the vegetation line
            # Calculate the distance from the vegetation line intersection point to the FS line
            distance_FS = intersection.distance(fs_line)

            if distance_FS > 500:
                alpha = (min(float(c['Db_max']),intersection_elevation) + float(c['FS_depth'])) / (beach_width+50) #This is a special case for when bathymetry is not available - alpha is calculated as the beach slope
            else:
                #Calclulate the slope from the vegetation line to the FS line
                alpha = (min(float(c['Db_max']),intersection_elevation) + float(c['FS_depth'])) / distance_FS
        
            # Append the morphological parameters to the lists
            s['Lveg'][index] = distance
            s['B'][index] = beach_width
            s['V'][index] = V
            s['h'][index] = h
            s['h_max'][index] = h_max
            s['alpha'][index] = alpha
            s['Db'][index] = intersection_elevation

    if p['process_plotting']:
        plot_sh_dir(s,sh_line)
            
    # Calculate morphological parameters for future scenarios
    if p['process_SLR']:
        #Calculate the distance between the road point and the historical coastline
        hist_distance_vegline(s,c,f)

        #Calculate the distance between the road points and the future coastline/vegline
        shoreline_change(s,c,f,p)


        #Iterate through road points to calculate morphological parameters for future scenarios
            # Iterate through the road computation points


        s['V_slr'] = [None]*len(s['x'])
        s['h_slr'] = [None]*len(s['x'])
        s['h_max_slr'] = [None]*len(s['x'])
        
        for index in range(len(s['x'])):
            if s['struc'][index]:
                s['V_slr'][index] = None
                s['h_slr'][index] = None
                s['h_max_slr'][index] = None
            else:
                road_point = Point(s['x'][index], s['y'][index])
                
                
                # Find the corresponding transect for the road point
                transect = transects.iloc[index].geometry
                elevations = transects.iloc[index].elevations
                
                if s['Lveg_slr'][index] > 0:
                    # Calculate the maximum height the sediment reaches above the intersection elevation with the future vegetation line
                    if s['Lveg_slr'][index] < s['Lveg'][index]:
                        # Filter the elevations between the road and future vegetation lines and subtract the dune toe height and slr from the elevations
                        start_index = int(c['road_width']) # Start half the road distance from the road point (Divided by two to get half the road width and multiplied by two because grid size i 0.5 m)
                        filtered_elevations = [(elev - (s['Db'][index]+float(c['SLR']))) for elev in elevations[start_index:int(s['Lveg_slr'][index] * 2)]]
                        if filtered_elevations:
                            # Calculate the maximum height of the filtered elevations
                            h_max_slr = max(filtered_elevations)
                            #Calculate the mean of the filtered elevations
                            h_slr = np.mean(filtered_elevations)
                            # Calculate the volume of sediment above the intersection elevation
                            V_slr = np.trapz(filtered_elevations, dx=0.5)
                    else:
                        V_slr = (s['h'][index]-float(c['SLR']))*s['Lveg_slr'][index]
                        h_slr = s['h'][index]-float(c['SLR'])
                        h_max_slr = s['h_max'][index]-float(c['SLR'])
                else :
                    V_slr = 0
                    h_max_slr = 0
                    h_slr = 0
                
                # Append the morphological parameters to the lists
                if V_slr < 0:
                    V_slr = 0
                    h_slr = 0
                    h_max_slr = 0
                
                s['V_slr'][index] = V_slr
                s['h_slr'][index] = h_slr
                s['h_max_slr'][index] = h_max_slr


    return s

def hist_distance_vegline(s,c,f):
    with shapefile.Reader(f['hist_CL_file']) as hist_cl:        
        # Extract the first shape (assuming there's only one polyline in the shapefile)
        hist_cl_shape = hist_cl.shape(0)
        hist_cl_points = hist_cl_shape.points

        # Create LineString objects from the points
        hist_cl_line = LineString(hist_cl_points)

        # Initialize a list to store the distance from the road to the historical vegetation line
        s['Lveg_hist'] = [None]*len(s['x'])

        # Iterate through the road computation points
        for index in range(len(s['x'])):
            road_point = Point(s['x'][index], s['y'][index])
            
            # Calculate the smallest distance to the historical vegetation line
            distance = road_point.distance(hist_cl_line)

            # Append the distance to the historical vegetation line to the list
            if s['struc'][index]:
                s['Lveg_hist'][index] = None
            else:
                s['Lveg_hist'][index] = distance

    return s

def bathyline_coords(s,c,f):
    with shapefile.Reader(f['bathyline_file']) as bathyline:
        # Extract the first shape (assuming there's only one polyline in the shapefile)
        bathyline_shape = bathyline.shape(0)
        bathyline_points = bathyline_shape.points

        # Create LineString objects from the points
        bathyline_line = LineString(bathyline_points)

        # Initialize lists to store the coordinates on the bathyline for wave extraction
        s['x_wave'] = [None]*len(s['x'])
        s['y_wave'] = [None]*len(s['x'])

        # Iterate through the road computation points
        for index in range(len(s['x'])):
            road_point = Point(s['x'][index], s['y'][index])
            
            # Find the nearest point on the bathyline
            nearest_bathy_point = nearest_points(road_point, bathyline_line)[1]

            # Append the distance to the historical vegetation line to the list
            s['x_wave'][index] = nearest_bathy_point.x
            s['y_wave'][index] = nearest_bathy_point.y

    if p['process_plotting']:
        plot_bathy_points(s,c,f)

    return s

# Plot the road and computation points
def plot_computation_points(s, c, f):
    # Open the shapefile
    sf = shapefile.Reader(f['road_file'])

    # Make a plot that shows the road and the computation points with DEM in the background and label the points with the z values
    with rasterio.open(f['DEM_file']) as dem:
        # Read the DEM data
        dem_data = dem.read(1)
        dem_extent = (dem.bounds.left, dem.bounds.right, dem.bounds.bottom, dem.bounds.top)
        
        # Make a plot that shows the road and the computation points with DEM in the background
        plt.figure()
        
        # Define the z value limits
        z_min = 0  # Set your desired minimum z value
        z_max = max(s['z']) + 1  # Set your desired maximum z value

        # Plot the DEM with z value limits
        im = plt.imshow(dem_data, extent=dem_extent, cmap='terrain', origin='upper', vmin=z_min, vmax=z_max)

        # Initialize variables to determine the extent of the road segment
        minx, miny, maxx, maxy = float('inf'), float('inf'), float('-inf'), float('-inf')

        # Iterate over all shapes in the shapefile
        for shape in sf.shapes():
            points = shape.points
            
            # Create a LineString object from the points
            line = LineString(points)

            # Plot the road
            plt.plot(*line.xy, color='black', linewidth=2)

            # Update the extent of the road segment
            line_minx, line_miny, line_maxx, line_maxy = line.bounds
            minx, miny = min(minx, line_minx), min(miny, line_miny)
            maxx, maxy = max(maxx, line_maxx), max(maxy, line_maxy)

        # Set aspect ratio to 1:1
        plt.gca().set_aspect('equal', adjustable='box')
        
        # Plot the computation points
        plt.scatter(s['x'], s['y'], color='red')
        
        # Label the points with the z values
        for x, y, z in zip(s['x'], s['y'], s['z']):
            plt.text(x, y, f'{z:.2f}', fontsize=8, ha='right', va='bottom')
        
        # Set labels and title
        plt.xlabel('X-coordinate [m]')
        plt.ylabel('Y-coordinate [m]')
        plt.title('Road and Computation Points with DEM Background')

        # Add a colorbar with the same scale as the DEM
        cbar = plt.colorbar(im)
        cbar.set_label('Elevation [m]')

        # Zoom in on the extent of the road segment
        buffer = 2 * float(c['spacing'])
        plt.xlim(minx - buffer, maxx + buffer)
        plt.ylim(miny - buffer, maxy + buffer)
        
        # Save the plot
        plt.savefig('out/road_and_computation_points_with_dem.png')


#plot transects and cl/veglines
def plot_transects(result_lines, cl_line, sh_line):
    # Plot the transects and the vegetation lines
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot the vegetation line
    gpd.GeoSeries([cl_line, sh_line]).plot(ax=ax, color='green', linewidth=2, label='Vegetation Line')
    
    # Plot the transects
    result_lines.plot(ax=ax, color='red', linewidth=1, label='Transects')
    
    # Set the aspect ratio to be equal
    ax.set_aspect('equal')
    
    # Set the labels and title
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Transects, Shoreline and Vegetation Lines')
    
    # Add a legend
    ax.legend()
    
    # Save the plot
    plt.savefig('out/transects_from_comp_points_shoreline.png')

    # Plot the elevation data along the transects, all transects in one plot
    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot the elevation data along the transects
    for i, elevations in enumerate(result_lines['elevations']):
        ax.plot(elevations, label=f'Transect {i+1}')


    # Set the labels and title
    ax.set_xlabel('Distance along Transect [m]')
    ax.set_ylabel('Elevation [m]')
    ax.set_title('Elevation Data along Transects')

    # Add a legend
    ax.legend()

    # Save the plot
    plt.savefig('out/elevation_data_along_transects.png')

def generate_transects(f):
        # Extract the vegetation line
    with shapefile.Reader(f['CL_file']) as cl, \
         shapefile.Reader(f['SH_file']) as sh, \
         shapefile.Reader(f['FS_file']) as fs:
        
        # Extract the first shape (assuming there's only one polyline in the shapefile)
        cl_shape = cl.shape(0)
        sh_shape = sh.shape(0)
        fs_shape = fs.shape(0)
        cl_points = cl_shape.points
        sh_points = sh_shape.points
        fs_points = fs_shape.points

        # Create LineString objects from the points
        cl_line = LineString(cl_points)
        sh_line = LineString(sh_points)
        fs_line = LineString(fs_points)

        crs_code = c['crs']
        
        # Create an empty GeoDataFrame to store the resulting lines
        result_lines = gpd.GeoDataFrame(columns=['geometry'],crs=crs_code)
        
    # Iterate through the road computation points
    for x, y in zip(s['x'], s['y']):
        road_point = Point(x, y)
            
        # Find the nearest point on the shoreline, sh_line
        nearest_shore_point = nearest_points(road_point, sh_line)[1]

        # Create a line from the road point to the nearest point
        optimized_line = LineString([road_point, nearest_shore_point])
            
        # Add the line to the result GeoDataFrame using pd.concat
        new_row = gpd.GeoDataFrame({'geometry': [optimized_line]}, crs=crs_code)
        result_lines = pd.concat([result_lines, new_row], ignore_index=True)

    # Extract elevation data from the DEM along the transects
    # Read the DEM file
    with rasterio.open(f['DEM_file']) as dem:
        dem_data = dem.read(1)
        dem_transform = dem.transform

        # Initialize a list to store the elevation data
        elevations = []

        # Iterate through the result lines
        for k, line in enumerate(result_lines['geometry'], start=1):

            # Calculate the total length of the line
            total_length = line.length

            # Generate distances at every 0.5 meters
            distances = np.arange(0, total_length, 0.5)

            # Interpolate points along the line at these distances
            points = [line.interpolate(distance) for distance in distances]

            # Create a new LineString with these points
            new_line = LineString(points)
            line = new_line

            # Extract the coordinates of the line
            coords = line.coords

            # Initialize a list to store the elevation data for the line
            line_elevations = []

            # Iterate through the coordinates
            for coord in coords:
                # Get the row and column indices corresponding to the x, y coordinates
                row, col = rasterio.transform.rowcol(dem_transform, coord[0], coord[1])
                # Read the elevation value at the row, col position
                z = dem_data[row, col]
                # Parse z result to float
                z = float(z)
                # Append the z value to the list
                line_elevations.append(z)

            # Append the list of z values to the elevations list
            elevations.append(line_elevations)


        # Add the elevation data to the result_lines GeoDataFrame
        result_lines['elevations'] = elevations


    if p['process_plotting']:
        plot_transects(result_lines, cl_line, sh_line)

    return result_lines, cl_line, sh_line, fs_line


def plot_bathy_points(s,c,f):
    # PLot bathypoints and bathylines and road and road points
    with shapefile.Reader(f['bathyline_file']) as bathyline, \
         shapefile.Reader(f['road_file']) as road:
        # Extract the first shape (assuming there's only one polyline in the shapefile)
        bathyline_shape = bathyline.shape(0)
        bathyline_points = bathyline_shape.points

        road_shape = road.shape(0)
        road_points = road_shape.points

        # Create LineString objects from the points
        bathyline_line = LineString(bathyline_points)
        road_line = LineString(road_points)

        # Make a plot with the road and bathypoints together with the bathyline and road line
        plt.figure()

        # Plot bathyline
        bathyline_x, bathyline_y = bathyline_line.xy
        plt.plot(bathyline_x, bathyline_y, label='Bathyline', color='blue')

        # Plot road line
        road_x, road_y = road_line.xy
        plt.plot(road_x, road_y, label='Road Line', color='green')

        # Plot s['x'], s['y']
        plt.scatter(s['x'], s['y'], label='s[x], s[y]', color='red')

        # Plot s['x_wave'], s['y_wave']
        plt.scatter(s['x_wave'], s['y_wave'], label='s[x_wave], s[y_wave]', color='orange')

        # Zoom in on the extent of the road segment
        minx, miny, maxx, maxy = road_line.bounds
        buffer = 2 * float(c['spacing'])
        plt.xlim(minx - buffer, maxx + buffer)
        plt.ylim(miny - buffer, maxy + buffer)

        # Add labels and legend
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.title('Plot of Bathyline, Road Line, and Points')
        plt.legend()

        # Save the plot
        plt.savefig('out/bathyline_wave_points.png')

# Function to calculate the orientation in degrees north
def calculate_orientation(line):
    dx = line.coords[1][0] - line.coords[0][0]
    dy = line.coords[1][1] - line.coords[0][1]
    angle = np.degrees(np.arctan2(dy, dx))
    orientation = (90 - angle) % 360
    return orientation        
        
def plot_sh_dir(s, sh_line):
        # Plot the roated lines and label them with the orientation
        fig, ax = plt.subplots(figsize=(10, 10))
        gpd.GeoSeries([sh_line]).plot(ax=ax, color='blue', linewidth=2)
        plt.scatter(s['x'], s['y'], color='red')


        # Label the orientation, check if orientation is None and handle it
        for x, y, sh_dir in zip(s['x'], s['y'], s['sh_dir']):
            if sh_dir is not None:
                plt.text(x, y, f'{sh_dir:.2f}', fontsize=8, ha='right', va='bottom')
            else:
                plt.text(x, y, 'N/A', fontsize=8, ha='right', va='bottom')

        # Zoom in on the road points
        minx, miny, maxx, maxy = min(s['x']), min(s['y']), max(s['x']), max(s['y'])
        buffer = 2 * float(c['spacing'])
        plt.xlim(minx - buffer, maxx + buffer)
        plt.ylim(miny - buffer, maxy + buffer)

        ax.set_aspect('equal')
        ax.set_xlabel('Easting')
        ax.set_ylabel('Northing')
        ax.set_title('Orientation of the Shoreline')
        plt.savefig('out/orientation_of_shoreline.png')

def check_struc(s,c,f):

    # open the shoreline shapefile
    with shapefile.Reader(f['SH_file']) as shoreline:
        # Extract the first shape (assuming there's only one polyline in the shapefile)
        shoreline_shape = shoreline.shape(0)
        shoreline_points = shoreline_shape.points

        # Create LineString objects from the points
        shoreline_line = LineString(shoreline_points)

    # open the non-erodible feature shapefile
    with shapefile.Reader(f['NE_file']) as ne:
        # Extract all shapes from the shapefile
        ne_shapes = ne.shapes()

        # Initialize a list to store the non-erodible features
        non_erodible_features = []

        # Iterate through the shapes
        for shape in ne_shapes:
            # Create LineString objects from the points
            non_erodible_features.append(LineString(shape.points))

    # Check if there is a non-erodible feature between the road and the shoreline (shortest distance from roadpoint)
    s['struc'] = [False]*len(s['x'])

    # Iterate through the road computation points
    for index in range(len(s['x'])):
        road_point = Point(s['x'][index], s['y'][index])
        
        # Find the nearest point on the shoreline
        nearest_shore_point = nearest_points(road_point, shoreline_line)[1]
        
        # Find the nearest non-erodible feature
        nearest_ne_feature = min(non_erodible_features, key=lambda feature: feature.distance(road_point))
        
        # Check if the non-erodible feature is between the road point and the shoreline
        if nearest_ne_feature.distance(road_point) < nearest_shore_point.distance(road_point):
            s['struc'][index] = True

    return s