# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 11:36:55 2023

@author: Caroline Hallin

This script contains functions to calculate the probability of flooding due to still water level (SWL) using the 
Generalized Extreme Value (GEV) distribution.

It reads SWL data from input files, fits the GEV distribution to the data, and calculates the probability of inundation 
dependning of the ground level in the different road points. In case of two SWL files, it interpolates the water levels along the road segment.

The impact of sea level rise (SLR) on the probability of inundation is accounted for by increasing the reference still water level.
"""

import numpy as np
from scipy.stats import genextreme
import matplotlib.pyplot as plt
import pandas as pd
import time as comp_time
import warnings
from numba import njit
from input import *
from dictionaries import *
from datetime import datetime, timedelta


# Fit GEV distribution to data
def GEV_SWL(s, c, f, p):
    
    if f['SWL_file_2'] is None:
        # Read swl input
        date, time, level = read_swl_input(f['SWL_file_1'])
    
        # Calculate the maximum level for each year
        year_max_level = year_max(date, time, level)

        # Fit the GEV distribution to the data
        params = genextreme.fit(year_max_level.Level)
        # Extract shape, location, and scale parameters
        shape, loc, scale = params

        if p['process_plotting']:

            # Generate a range of values for the x-axis
            x = np.linspace(min(year_max_level.Level), max(year_max_level.Level), num=100)
            # Calculate the probability density function (PDF) using the fitted parameters
            pdf = genextreme.pdf(x, shape, loc, scale)
                    # Plot the histogram of the data and the fitted GEV PDF
            plt.figure()
            plt.hist(year_max_level.Level, bins=20, density=True, alpha=0.6, label="Data Histogram")
            plt.plot(x, pdf, 'r', label="Fitted GEV PDF")
            plt.xlabel("Yearly maximum water level [m rel MSL]")
            plt.ylabel("Probability Density")
            plt.title("Fitting GEV Distribution")
            plt.legend()

            # Save the plot to the 'out' folder
            plt.savefig('out/gev_distribution.png')
            plt.close()

            # Calculate the return levels for different return periods
            return_periods_eq = np.linspace(1, 100, num=100)
            return_levels = genextreme.ppf(1 - 1 / np.array(return_periods_eq), shape, loc, scale)
            # Observed data
            levels = year_max_level.Level

            # Step 1: Calculate the ranks
            ranks = np.argsort(np.argsort(-levels)) + 1

            # Step 2: Calculate the return periods
            N = len(levels)
            return_periods = (N + 1) / ranks


            # Plot the return levels
            plt.figure()
            #Plot the observed data
            plt.plot(return_periods, levels, 'o')
            plt.plot(return_periods_eq, return_levels, 'b-')
            plt.xlabel("Return Period (years)")
            plt.ylabel("Return Level [m rel MSL]")
            plt.title("GEV SWL")
            plt.xlim([0, 100])
            # Save the plot to the 'out' folder
            plt.savefig('out/return_levels_SWL.png')
            plt.close()

        s['p_swl'] = prob_swl_one(s, shape, loc, scale, float(c['MSL']))

        # Calculate the 100-year return level
        s['swl_100'] = [genextreme.ppf(1 - 1 / 100, shape, loc, scale)]* len(s['ID'])

        if p['process_SLR']:
            # Calculate the probability of flooding due to SWL after SLR
            s['p_swl_slr'] = prob_swl_one(s, shape, loc, scale, float(c['MSL']) + float(c['SLR']))

    else:
        s['p_swl'] = [None] * len(s['ID'])
        s['p_swl_slr'] = [None] * len(s['ID'])
        s['swl_100'] = [None] * len(s['ID'])
        

        date1_num, time1_num, level1, date2_num, time2_num, level2, reference_date = read_two_SWL_files(f['SWL_file_1'], f['SWL_file_2'])

        
        # Interpolate the water levels for the same dates
        for i in range(len(s['ID'])):

            date, time, level = [], [], []
            x = s['x'][i]
            y = s['y'][i]
            x1 = float(c['SWL_1_x'])
            y1 = float(c['SWL_1_y'])
            x2 = float(c['SWL_2_x'])
            y2 = float(c['SWL_2_y'])
            

            # Run the interpolation
            date_num, time_num, level = interpolate_water_levels(date1_num, time1_num, level1, date2_num, time2_num, level2, x, y, x1, y1, x2, y2)
            
            # Convert numerical dates back to original format
            date = [(reference_date + timedelta(days=int(d))).strftime('%Y-%m-%d') for d in date_num]
            
            # Convert numerical times back to original format
            time = [str(timedelta(seconds=int(t))) for t in time_num]

            # Calculate the maximum level for each year
            year_max_level = year_max(date, time, level)

            # Fit the GEV distribution to the data
            params = genextreme.fit(year_max_level.Level)
            # Extract shape, location, and scale parameters
            shape, loc, scale = params

            # Calculate the 100-year return level
            s['swl_100'][i] = genextreme.ppf(1 - 1 / 100, shape, loc, scale)

            s['p_swl'][i] = prob_swl(s, shape, loc, scale, float(c['MSL']),i)

            if p['process_plotting']:

                # Generate a range of values for the x-axis
                x = np.linspace(min(year_max_level.Level), max(year_max_level.Level), num=100)
                # Calculate the probability density function (PDF) using the fitted parameters
                pdf = genextreme.pdf(x, shape, loc, scale)
                        # Plot the histogram of the data and the fitted GEV PDF
                plt.figure()
                plt.hist(year_max_level.Level, bins=20, density=True, alpha=0.6, label="Data Histogram")
                plt.plot(x, pdf, 'r', label="Fitted GEV PDF")
                plt.xlabel("Yearly maximum water level [m rel MSL]")
                plt.ylabel("Probability Density")
                plt.title("Fitting GEV Distribution")
                plt.legend()

                # Save the plot to the 'out' folder
                plt.savefig(f'out/gev_distribution{s["ID"][i]}.png')
                plt.close()

                # Calculate the return levels for different return periods
                return_periods_eq = np.linspace(0, 100, num=100)
                return_levels = genextreme.ppf(1 - 1 / np.array(return_periods_eq), shape, loc, scale)
                # Observed data
                levels = year_max_level.Level

                # Step 1: Calculate the ranks
                ranks = np.argsort(np.argsort(-levels)) + 1

                # Step 2: Calculate the return periods
                N = len(levels)
                return_periods = (N + 1) / ranks


                # Plot the return levels
                plt.figure(figsize=(3.5, 3))
                #Plot the observed data
                plt.plot(return_periods, levels, 'bo')
                plt.plot(return_periods_eq, return_levels, 'k-')
                plt.plot(return_periods_eq, return_levels+float(c['SLR']), 'k--')
                plt.legend(['Observations', rf"Fitted GEV ($\mu = {loc:.2f}$, $\sigma = {scale:.2f}$, $\xi = {shape:.2f}$)", 'Fitted GEV + SLR'], loc='upper center', bbox_to_anchor=(0.4, -0.20), ncol=1)
                plt.xlabel("Return period (years)")
                plt.ylabel("Return level (m rel MSL)")
                plt.title("Still water level")
                plt.xlim([0, 100])
                # Save the plot to the 'out' folder
                plt.savefig(f'out/return_levels{s["ID"][i]}.png',dpi=300, bbox_inches='tight')
                plt.close()


            if p['process_SLR']:
                # Calculate the probability of flooding due to SWL after SLR
                s['p_swl_slr'][i] = prob_swl(s, shape, loc, scale, float(c['MSL']) + float(c['SLR']),i)     
        
    return s

# Read SWL input   
def read_swl_input(SWL_file):
    # Open the SWL file and read the time and level data
    with open(SWL_file, 'r') as file:
        data = file.readlines()
        date = []
        time = []
        level = []
        for line in data:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            date.append(parts[0])
            time.append(parts[1])
            level.append(float(parts[2]))
    return date, time, level

# Calculate the maximum level for each year
def year_max(date, time, level):
    # Create a DataFrame with the date and level data
    df = pd.DataFrame({'Date': date, 'Time': time, 'Level': level})

    # Convert the 'Date' column to datetime objects
    df['Date'] = pd.to_datetime(df['Date'])
    
    # Shift year boundary to July
    df['Year'] = df['Date'].apply(lambda x: x.year if x.month < 7 else x.year + 1)

    # Calculate the average water level for each year
    year_avg = df.groupby('Year')['Level'].mean().reset_index()
    
    # Group by 'Year' and find the maximum 'Level' for each year
    yearly_max = df.groupby('Year')['Level'].max().reset_index()

    # Subtract the average water level from the maximum level for each year
    yearly_max['Level'] = yearly_max['Level'] - year_avg['Level']

    # Return the maximum level for each year
    return yearly_max

# Function to calculate the probability of flooding due to SWL
def prob_swl(s, shape, loc, scale, ref,i):
    # Calculate the return period for inundation of the computation point

    # Calculate the return level for each computation point
    return_level = s['z'][i] - ref

    # Calculate the return period for the computed return level
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        return_period = 1 / (1 - genextreme.cdf(return_level, shape, loc, scale))
    
    return return_period

# Function to calculate the probability of flooding due to SWL
def prob_swl_one(s, shape, loc, scale, ref):
    # Calculate the return period for inundation of the computation point

    # Calculate the return level for each computation point
    return_level = [z - ref for z in s['z']]

    # Calculate the return period for the computed return level
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        return_period = 1 / (1 - genextreme.cdf(return_level, shape, loc, scale))
    
    return return_period


def interpolate_water_levels(date1_num, time1_num, level1, date2_num, time2_num, level2, x, y, x1, y1, x2, y2):
    # Combine date and time into a single numerical value for interpolation
    datetime1 = date1_num + time1_num / (24 * 3600)
    datetime2 = date2_num + time2_num / (24 * 3600)
    
    # Find overlapping time steps
    start_time = max(datetime1[0], datetime2[0])
    end_time = min(datetime1[-1], datetime2[-1])
    
    # Create a mask for the overlapping time steps
    mask1 = (datetime1 >= start_time) & (datetime1 <= end_time)
    mask2 = (datetime2 >= start_time) & (datetime2 <= end_time)
    
    # Filter the arrays based on the mask
    datetime1_filtered = datetime1[mask1]
    datetime2_filtered = datetime2[mask2]
    level1_filtered = level1[mask1]
    level2_filtered = level2[mask2]
    date1_num_filtered = date1_num[mask1]
    time1_num_filtered = time1_num[mask1]
    
    # Find common time steps
    common_times, idx1, idx2 = np.intersect1d(datetime1_filtered, datetime2_filtered, return_indices=True)
    
    # Initialize arrays to store interpolated results
    interpolated_dates = date1_num_filtered[idx1]
    interpolated_times = time1_num_filtered[idx1]
    interpolated_levels = np.empty(len(common_times), dtype=level1.dtype)
    
    # Calculate the distance between the two points
    dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    if dist == 0:
        # Handle the case where the distance is zero
        interpolated_levels = (level1_filtered[idx1] + level2_filtered[idx2]) / 2
    else:
        # Calculate the weights for interpolation
        weight1 = np.sqrt((x2 - x)**2 + (y2 - y)**2) / dist
        weight2 = np.sqrt((x1 - x)**2 + (y1 - y)**2) / dist
        
        # Interpolate the water levels
        interpolated_levels = (weight1 * level1_filtered[idx1] + weight2 * level2_filtered[idx2]) / (weight1 + weight2)
    
    return interpolated_dates, interpolated_times, interpolated_levels

def read_two_SWL_files(SWL_file_1, SWL_file_2):
    # Read the SWL input files
    date1, time1, level1 = read_swl_input(SWL_file_1)
    date2, time2, level2 = read_swl_input(SWL_file_2)

    # Define a reference date
    reference_date = datetime(1900, 1, 1)

    # Convert date strings to datetime objects and then to numerical format
    date1_num = np.array([(datetime.strptime(d, '%Y-%m-%d') - reference_date).days for d in date1], dtype=np.float64)
    date2_num = np.array([(datetime.strptime(d, '%Y-%m-%d') - reference_date).days for d in date2], dtype=np.float64)

    # Convert time strings to numerical format (seconds since midnight)
    time1_num = np.array([int(t.split(':')[0]) * 3600 + int(t.split(':')[1]) * 60 + int(t.split(':')[2]) for t in time1], dtype=np.float64)
    time2_num = np.array([int(t.split(':')[0]) * 3600 + int(t.split(':')[1]) * 60 + int(t.split(':')[2]) for t in time2], dtype=np.float64)

    # Ensure level inputs are NumPy arrays with appropriate data types
    level1 = np.array(level1, dtype=np.float64)
    level2 = np.array(level2, dtype=np.float64)

    return date1_num, time1_num, level1, date2_num, time2_num, level2, reference_date