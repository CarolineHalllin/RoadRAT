# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 09:17:03 2023

@author: Caroline Hallin

Script to calculate the long-term shoreline change due to sea level rise (SLR) and projected observed long-term erosion/accretion.

"""
import numpy as np
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point
import xarray as xr
import pandas as pd
import pyproj
import glob
import re
import time
from dictionaries import *
from runup_erosion import rev_shoal

# Calculate the shoreline change due to SLR and longterm erosion/accretion
def shoreline_change(s,c,f,p):
    s['Lveg_slr'] = [None]*len(s['ID'])
    #Calculate depth of closure
    if p['process_extract_waves']:
        closure_depth_extract_waves(s,c,f)
    else:
        closure_depth_wave_file(s,c,f)

    #Calculate background erosion/accretion in m/year
    sl_change_rate = np.zeros(len(s['ID']))
    for i in range(len(s['ID'])):
        if s['struc'][i]:
            s['Lveg_slr'][i] = None
        else:
            sl_change_rate[i]=(s['Lveg_hist'][i]-s['Lveg'][i])/(int(c['start_yr'])-int(c['CL_yr']))

    # Calculate length of Dean profile
    dean(s,c,f)
            
    #Calculate shoreline change due to SLR and longterm erosion/accretion

    for i in range(len(s['ID'])):
        if s['struc'][i]:
            s['Lveg_slr'][i] = None
        else:
            s['Lveg_slr'][i] = s['Lveg'][i] - float(c['SLR']) * s['LDean'][i] / s['Dc'][i] - sl_change_rate[i]*(int(c['end_yr'])-int(c['start_yr']))
            
    return s


def closure_depth_extract_waves(s, c, f):

    s['Dc'] = [None] * len(s['ID'])
    wave_dir = r"{}".format(f['wave_files_dir'])

    # Convert wave input coordinates to the same coordinate system as the wave data
    wave_points = [Point(x, y) for x, y in zip(s['x_wave'], s['y_wave'])]

    crs_point = extract_numbers(c['crs'])
    crs_wave = extract_numbers(c['crs_wave'])

    points_df = gpd.GeoDataFrame(geometry=wave_points, crs=int(crs_point[0]))
    points_df_wgs84 = points_df.to_crs(epsg=int(crs_wave[0]))

    points = [list(row['geometry'].coords)[0] for index, row in points_df_wgs84.iterrows()]

    # Find the wave data file index corresponding to the closest point of the wave extraction point
    first_file = glob.glob(wave_dir + r"\*.nc")[0]
    ds = xr.open_dataset(first_file, engine='netcdf4')

    geod = pyproj.Geod(ellps='WGS84')
    coord = pd.DataFrame({'x': ds.longitude[:], 'y': ds.latitude[:]})

    # Find the unique indices
    closest_indices = []
    unique_indices = []

    for row in points:
        lat, lon = row[1], row[0]
        distances = [geod.inv(coord['x'].iloc[i], coord['y'].iloc[i], lon, lat)[2] for i in range(len(coord))]
        closest_index = int(np.argmin(distances))
        closest_indices.append(closest_index)
        if closest_index not in unique_indices:
            unique_indices.append(closest_index)

    # Initialize array to store depth of closure
    wave_files = glob.glob(wave_dir + r"\*.nc")
    dc = np.zeros((len(unique_indices), len(wave_files)))

    # Read wave data
    for k, file in enumerate(wave_files):
        ds = xr.open_dataset(file, engine='netcdf4')

        for j, index in enumerate(unique_indices):
            model_data_temp = pd.DataFrame({'Date': ds.time[264:], 'depth': ds.botl[264:, index],
                                            'Hm0': ds.hs[264:, index], 'tp': ds.tps[264:, index]})

            model_data_temp.set_index(['Date'], inplace=True)

            # Reverse shoaling to calculate the deep water wave height
            Hin = model_data_temp['Hm0'].to_numpy()
            Tin = model_data_temp['tp'].to_numpy()
            din = model_data_temp['depth'].to_numpy()
            theta = np.zeros(len(Hin))

            model_data_temp['Hm0'] = rev_shoal(Tin, din, Hin, theta)

            # Hallermeier
            temp = model_data_temp.sort_values(by=['Hm0'], ascending=False)
            H = temp['Hm0'].iloc[11]  # Adjust if wave data is not hourly
            Tp = temp['tp'].iloc[11]
            dc_value = 2.28 * H - 68.5 * H**2 / (9.81 * Tp**2)
            dc[j, k] = dc_value

    dc_ave = dc.mean(axis=1)

    # Initialize S['Dc'] and populate with the average depth of closure
    for i in range(len(s['ID'])):
        if s['struc'][i]:
            s['Dc'][i] = None
        else:
            s['Dc'][i] = dc_ave[unique_indices.index(closest_indices[i])]


    return s


def closure_depth_wave_file(s,c,f):
    #Read wave file
    wave_file = f['wave_file']
    wave_data = pd.read_csv(wave_file, sep=',')
    # Convert 'datetime' column to datetime object
    wave_data['datetime'] = pd.to_datetime(wave_data['datetime'])

    # Set 'datetime' as the index
    wave_data.set_index('datetime', inplace=True)

    # Determine the range of years
    min_year = wave_data.index.year.min()
    max_year = wave_data.index.year.max()

    # Initialize arrays to store results
    num_years = max_year - min_year + 1
    dc = np.zeros(num_years)
    
    # Group by year
    grouped = wave_data.groupby(wave_data.index.year)

    
    for year, group in grouped:
        # Sort wave heights in descending order
        sorted_data = group.sort_values(by='significant_wave_height_m', ascending=False)

        # Calculate the number of samples corresponding to 12 hours
        num_samples_12h = int(12 * 60 / ((group.index[1] - group.index[0]).seconds / 60))

        # Get the wave height and period at the 12-hour mark
        H_12h = sorted_data['significant_wave_height_m'].iloc[num_samples_12h - 1]
        T_12h = sorted_data['peak_wave_period_s'].iloc[num_samples_12h - 1]

        # Calculate the depth of closure
        dc[year - min_year]= 2.28 * H_12h - 68.5 * H_12h ** 2 / (9.81 * T_12h ** 2)

    # Initialize S['Dc'] and populate with the average depth of closure
    dc_ave=dc.mean()
    s['Dc'] = [dc_ave]*len(s['ID'])

    for i in range(len(s['ID'])):

        if s['struc'][i]:
            s['Dc'][i] = None

    return s


def dean(s,c,f):
    if f['grain_size_file']:
        read_grain_size(s,f)
    else:
        s['d50'] = [c['d50_default']]*len(s['ID'])

    spec_grav = 2.65 # Specific gravity of sediment
    cd = 10.9 # Drag coefficient
    g = 9.81 # Acceleration due to gravity

    s['LDean'] = [None]*len(s['ID'])

    for i in range(len(s['ID'])):
        if s['struc'][i]:
            s['LDean'][i] = None
        else:
            D_star=(g*(spec_grav-1)/c['viscosity']**2)**(1/3)*s['d50'][i]
            ws = c['viscosity']/s['d50'][i]*(np.sqrt(10.36**2+1.049*D_star**3)-10.36)
            A = 2.25*(ws**2/g)**(1/3)
            s['LDean'][i] = (s['Dc'][i]/A)**(3/2)

    return s


def read_grain_size(s,f):
    #add script to read grain size file
    return s

def extract_numbers(input_string):
    # Find all sequences of digits in the input string
    numbers = re.findall(r'\d+', input_string)
    # Convert the found sequences to integers
    return [int(num) for num in numbers]