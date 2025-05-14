# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 11:36:55 2023

@author: Caroline Hallin

Script to calculate the probability of runup reaching the road and the probability of storm erosion
reaching the road in present and future climate scenarios.

"""

import glob
import numpy as np
import pandas as pd
import pyproj
import xarray as xr
import geopandas as gpd
import re
import matplotlib.pyplot as plt
import time as time_compute
import logging
import warnings
from numba import jit
from shapely.geometry import Point
from scipy.stats import genextreme
from scipy.stats import genpareto
from SWL import read_swl_input, interpolate_water_levels, read_two_SWL_files
from datetime import datetime, timedelta
from scipy.optimize import curve_fit
from scipy.stats import rankdata


# Calculate the probability of runup reaching the road
def prob_runup_erosion(s,c,f,p):
    
    if p['process_runup'] or p['process_erosion']:
        s['p_r'] = np.zeros(len(s['ID']))
        s['Db_r'] = [None] * len(s['ID'])
        s['V_r'] = [None] * len(s['ID'])
        s['r_100'] = [None] * len(s['ID'])
        if p['process_SLR']:
            s['p_r_slr'] = np.zeros(len(s['ID']))
            s['V_r_slr'] = [None] * len(s['ID'])
    if p['process_erosion']:
        s['p_er'] = np.zeros(len(s['ID']))
        s['er_100'] = [None] * len(s['ID'])
        if p['process_SLR']:
            s['p_er_slr'] = np.zeros(len(s['ID']))
            s['er_slr_100'] = [None] * len(s['ID'])

    # Define c_impact erosion coefficient for each computation point
    if p['process_erosion']:
        s['c_impact'] = np.zeros(len(s['ID']))
        if f['erodibility_file'] is None:
                s['c_impact'][:]=c['c_impact_default']
        
            #Add code here to read and extract C_impact coefficient from shaefile- actually the reading should take place directly under 
    
    #constants
    g=9.81

    if f['SWL_file_2'] is None:
        # Read water level data from file
        date, time, level = read_swl_input(f['SWL_file_1'])

        # Subtract the yearly average water level from the actual level
        df_swl = year_ave(date, time, level)
    else:
        # Read the SWL input files
        date1_num, time1_num, level1, date2_num, time2_num, level2, reference_date = read_two_SWL_files(f['SWL_file_1'], f['SWL_file_2'])
    
    # Loop over each computation point and calculate the probability of runup reaching the road
    
    #Check if there is a non-erodible layer between the road and the shoreline (struc = 'y')
    for i in range(len(s['ID'])):
        logging.info(f'Processing runup/erosion point ID: {s["ID"][i]}')

        if s['struc'][i]: #Don't calculate runup and erosion if there is a non-erodible layer between the road and the shoreline
            if p['process_runup'] or p['process_erosion']:
                s['p_r'][i] = None
                s['Db_r'][i] = None
                s['V_r'][i] = None
                s['r_100'][i] = None
                if p['process_SLR']:
                    s['p_r_slr'][i] = None
                    s['V_r_slr'][i] = None
            if p['process_erosion']:
                s['p_er'][i] = None
                s['er_100'][i] = None
                if p['process_SLR']:
                    s['p_er_slr'][i] = None
                    s['er_slr_100'][i] = None
        else:


            if p['process_extract_waves']:
                df_waves = extract_wave_data(s, c, f, i) #Read data from several netCDF files
            else:
                df_waves = read_wave_data(s, c, f, i) #Read data from a single csv file

            if f['SWL_file_2'] is not None:
            
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

                df_swl = year_ave(date, time, level)
            
            
            # Find wave and water level data that match the same date and time
            # Merge the wave and water level data on the 'Date' column
            df = pd.merge(df_waves, df_swl, on='Date', how='inner')

            # Determine time step in df
            dt = df['Date'].diff().dt.total_seconds()
            
            #Calculate wave height and angle (in radians) at breaking
            H = df['Hm0'].to_numpy()
            T= df['tp'].to_numpy()
            d=df['depth'].to_numpy()
            theta = df['theta'].to_numpy()
            dt=dt.to_numpy()

            Hb, dirb, db = larson_direct_method(H,T,d,theta,s['sh_dir'][i])

            df['TWL'] = np.zeros(len(df))
            
            # Reverse shoaling to find the equivalent wave height at deep water

            Hprim_onshore = rev_shoal(T,db,Hb,dirb)

            SWL=df['Level'].to_numpy()

            if p['process_runup'] or p['process_erosion']:

                # Calculate runup
                
                TWL = stockdon(Hprim_onshore, T, s['alpha'][i],SWL)

                df['TWL'] = TWL

                # GEV analysis on the runup data
                #Find year max
                year_max_TWL = year_max_total_WL(df)

                #Fit a GEV distribution to the data
                params = genextreme.fit(year_max_TWL)
                #Extract shape, location, and scale parameters
                shape, loc, scale = params

                #Calculate the 100-year total water level
                s['r_100'][i] = genextreme.ppf(1 - 1 / 100, shape, loc, scale)

                #Define dune toe elevation based on runup level with 10 year return period
                x = 10
                s['Db_r'][i] = genextreme.ppf(1 - 1 / x, shape, loc, scale) + float(c['MSL'])

                # Recalculate volume of sediment between the road and the dune toe
                if s['h'][i] > 0:
                    s['V_r'][i] = s['V'][i]+s['V'][i]/s['h'][i]*(s['Db'][i]-min(float(c['Db_max']),s['Db_r'][i]))
                else:
                    s['V_r'][i] = 0

                if s['V_r'][i] < 0:
                    s['V_r'][i] = 0

                if p['process_SLR']:
                     if s['h'][i] + (s['Db'][i]-s['Db_r'][i]) - float(c['SLR']) > 0:
                        s['V_r_slr'][i] = s['V_slr'][i]+s['V_slr'][i]/(s['h'][i]-float(c['SLR']))*(s['Db'][i]-min(float(c['Db_max']),s['Db_r'][i]))
                     else:
                         s['V_r_slr'][i] = 0

                     if s['V_r_slr'][i] < 0:
                         s['V_r_slr'][i] = 0
                

                if p['process_plotting']:
                    # Plot the GEV distribution
                    plot_GEV_TWL(year_max_TWL, shape, loc, scale, i,s,c)

                # Define the return level
                #Calculate the return level for each computation point
                MSL=float(c['MSL'])
                return_level = max((s['h_max'][i]+s['Db'][i]-MSL),(s['z'][i]-MSL))

                #Calculate the return period for the computed return level
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=RuntimeWarning)
                    return_period = 1/(1 - genextreme.cdf(return_level, shape, loc, scale))

                s['p_r'][i] = return_period

                # Calculate the probability of runup reaching the road in the future
                if p['process_SLR']:
                    SLR= float(c['SLR'])
                    # Calculate the return level for each computation point in the future
                    return_level_slr = max((s['h_max_slr'][i]+s['Db'][i]-MSL),(s['z'][i]-MSL-SLR))

                    # Calculate the return period for the computed return level
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", category=RuntimeWarning)
                        return_period_slr = 1/(1 - genextreme.cdf(return_level_slr, shape, loc, scale))

                    s['p_r_slr'][i] = return_period_slr

            # CaLculate storm erosion based on the Larson impact model
            if p['process_erosion']:
                zb = s['Db_r'][i] - float(c['MSL']) 
                ave_height = s['h'][i]
                V=s['V'][i]
                if ave_height < 0:
                    ave_height = 0
                
                c_imp=s['c_impact'][i]
                dist = s['B'][i]
                slope = (s['Db'][i]-float(c['MSL']))/s['B'][i]
                L_veg = s['Lveg'][i]
                                
                storm_erosion_dt = storm_erosion(zb,Hprim_onshore,T,SWL,c_imp,dt,dist,ave_height,slope,L_veg,V)



                if p['process_calibration']:
                    #Calculate the sum of erosion during calibration period
                    sum_erosion = sum(storm_erosion_dt)
                    #Print the point ID and sum of erosion during calibration period
                    print(f'Point ID: {s["ID"][i]}, Sum of erosion during calibration period: {sum_erosion} m3')
                    continue

                df['erosion'] = storm_erosion_dt

            # Calculated the annual erosion volume
            yearly_erosion = sum_erosion_per_year(df)
            threshold = 0.1
            excesses = yearly_erosion[yearly_erosion > threshold]

            # Check if excesses is not empty
            if len(excesses) == 0:
                if s['V'][i]==0: #s['V_r'][i] == 0:
                    s['p_er'][i] = 0
                    s['er_100'][i] = 0
                else:
                    s['p_er'][i] = np.inf
                    s['er_100'][i] = 0
            else:
                if p['process_hallermeier_rhodes']:
                    # Fit the Hallermeier and Rhodes equation to the storm erosion data
                    # Rank the data
                    excesses = np.sort(excesses)[::-1]
                    # Assign a rank to each data point
                    ranks = rankdata(-excesses, method='max')
                    # Adjust the ranks to the lentgh of the time sereíes
                    ranks = len(yearly_erosion) - ranks + 1
                    # Calculate the return period for each data point
                    return_periods = 1 / (1 - ranks / (len(yearly_erosion) + 1))

                    # Fit the Hallermeier and Rhodes equation to the data V=c*T^0.4 -> y=a*x^0.4
                    # Define model function 
                    def model(x, a):
                        return a * x ** 0.4
                    
                    # Fit the model to the data
                    params, _ = curve_fit(model, return_periods, excesses)
                    # Extract the parameter
                    a = params[0]

                

                    if p['process_plotting']:
                        # Plot the Hallermeier and Rhodes equation
                        return_periods_eq = np.linspace(0, 100, num=100)
                        plt.figure(figsize=(3.5, 3))
                        plt.plot(return_periods, excesses, 'ro', label="Simulated storm erosion")
                        plt.plot(return_periods_eq, model(return_periods_eq, a), 'k', label=rf"Hallermeier and Rhodes Eq ($c = {a:.2f}$)")
                        plt.xlabel("Return period (years)")
                        plt.xlim([0, 100])
                        plt.ylabel("Return volume (m$^3$/m)")
                        plt.title("Annual storm erosion")
                        plt.legend(loc='upper center', bbox_to_anchor=(0.4, -0.20), ncol=1)
                        # Save the plot to the 'out' folder
                        plt.savefig(f'out/hallermeier_rhodes_{s["ID"][i]}.png', dpi=300, bbox_inches='tight')
                        plt.close()

                    # Calculate the return period for the erosion volume exceeding the volume of sediment between the veg line and the road
                    return_level_erosion = s['V'][i]
                    
                    if return_level_erosion == 0:
                        s['p_er'][i] = 0
                    else:
                        s['p_er'][i] = (return_level_erosion / a) ** (1 / 0.4)

                    # Calculate the 100-year erosion volume
                    s['er_100'][i] = model(100, a)


                #Do the same for the future SLR
                if p['process_SLR']:
                    zb = s['Db_r'][i] - float(c['MSL']) 
                    ave_height = s['h_slr'][i]
                    V = s['V_slr'][i]
                    if ave_height < 0:
                        ave_height = 0
                    
                    c_imp=s['c_impact'][i]
                    dist = s['B'][i]
                    slope = (s['Db'][i]-float(c['MSL']))/s['B'][i]
                    L_veg = s['Lveg_slr'][i]

                    if L_veg < 0:
                        storm_erosion_dt_slr = np.zeros(len(df))
                    else:
                        storm_erosion_dt_slr = storm_erosion(zb,Hprim_onshore,T,SWL,c_imp,dt,dist,ave_height,slope,L_veg,V)

                   

                    df['erosion'] = storm_erosion_dt_slr

                # Calculated the annual erosion volume
                yearly_erosion = sum_erosion_per_year(df)
                threshold = 0.1
                excesses = yearly_erosion[yearly_erosion > threshold]


                # Check if excesses is not empty
                if len(excesses) == 0:
                    if s['V_slr'][i] == 0:
                        s['p_er_slr'][i] = 0
                        s['er_slr_100'][i] = 0
                    else:
                        s['p_er_slr'][i] = np.inf
                        s['er_slr_100'][i] = 0
                        
                else:

                    if p['process_hallermeier_rhodes']:
                        # Fit the Hallermeier and Rhodes equation to the storm erosion data
                        # Rank the data
                        excesses = np.sort(excesses)[::-1]
                        # Assign a rank to each data point
                        ranks = rankdata(-excesses, method='max')
                        # Adjust the ranks to the lentgh of the time sereíes
                        ranks = len(yearly_erosion) - ranks + 1
                        # Calculate the return period for each data point
                        return_periods = 1 / (1 - ranks / (len(yearly_erosion) + 1))
                        #print (excesses)

                        # Fit the Hallermeier and Rhodes equation to the data V=c*T^0.4 -> y=a*x^0.4
                        # Define model function 
                        def model(x, a):
                            return a * x ** 0.4
                        
                        # Fit the model to the data
                        params, _ = curve_fit(model, return_periods, excesses)
                        # Extract the parameter
                        a = params[0]

                    

                        if p['process_plotting']:
                            # Plot the Hallermeier and Rhodes equation
                            return_periods_eq = np.linspace(0, 100, num=100)
                            plt.figure(figsize=(3.5, 3))
                            plt.plot(return_periods, excesses, 'ro', label="Simulated storm erosion")
                            plt.plot(return_periods_eq, model(return_periods_eq, a), 'k', label=rf"Hallermeier and Rhodes Eq ($c = {a:.2f}$)")
                            plt.xlabel("Return period (years)")
                            plt.xlim([0, 100])
                            plt.ylabel("Return volume (m$^3$/m)")
                            plt.title("Annual storm erosion")
                            plt.legend(loc='upper center', bbox_to_anchor=(0.4, -0.20), ncol=1)
                            # Save the plot to the 'out' folder
                            plt.savefig(f'out/hallermeier_rhodes_SLR_{s["ID"][i]}.png', dpi=300, bbox_inches='tight')
                            plt.close()

                        # Calculate the return period for the erosion volume exceeding the volume of sediment between the veg line and the road
                        return_level_erosion_slr = s['V_slr'][i]
                        
                        if return_level_erosion_slr == 0:
                            s['p_er_slr'][i] = 0
                        else:
                            s['p_er_slr'][i] = (return_level_erosion_slr / a) ** (1 / 0.4)

                        # Calculate the 100-year erosion volume
                        s['er_slr_100'][i] = model(100, a)
                    

                    else:


                        # Fit a GPD to the storm erosion data
                        params = genpareto.fit(excesses)
                        # Extract the shape, location, and scale parameters
                        shape, loc, scale = params


                        if p['process_plotting']:
                            # Plot the GPD distribution
                            
                            plot_GPD_erosion(yearly_erosion, excesses, shape, loc, scale, i)

                        # Calculate the probabilty of the erosion volume exceeding the volume of sediment between the veg line and the road
                        num_exceedances = len(excesses)
                        total_data_points = len(yearly_erosion)
                        exceedance_prob = num_exceedances / total_data_points

                        #Calculate the return period for a specific return level
                        return_level_erosion = s['V'][i]
                        return_period_erosion = 1/(exceedance_prob*(1 - genpareto.cdf(return_level_erosion, shape, loc, scale)))

                        s['p_er'][i] = return_period_erosion

                        if p['process_SLR']:
                            # Calculate the return period for the erosion volume in the future
                            return_level_erosion_slr = s['V_slr'][i]
                            return_period_erosion_slr = 1/(exceedance_prob*(1 - genpareto.cdf(return_level_erosion_slr, shape, loc, scale)))

                            s['p_er_slr'][i] = return_period_erosion_slr

    return s


# Calculate the average water level for each year
def year_ave(date, time, level):
    # Create a DataFrame with the date and level data
    df = pd.DataFrame({'Date': date, 'Time': time, 'Level': level})

    # Parse the 'Date' and 'Time' columns to a single datetime column
    df['Date'] = pd.to_datetime(df['Date'] + ' ' + df['Time'])
    df['Date'] = pd.to_datetime(df['Date'])  # Convert 'Date' column to datetime
    df.set_index('Date', inplace=True)       # Set 'Date' column as the index

    # Drop the 'Time' column
    df.drop(columns='Time', inplace=True)

    # Calculate the average water level for each year
    year_avg = df.groupby(df.index.year)['Level'].mean()

    # Loop through each year and subtract the average water level
    for year in year_avg.index:
        df.loc[df.index.year == year, 'Level'] -= year_avg[year]

    # Return the data frame with the adjusted water levels
    return df

def extract_wave_data(s, c, f, i):
    # Read wave data
    wave_dir = r"{}".format(f['wave_files_dir'])

    # Convert wave input coordinates to the same coordinate system as the wave data
    wave_point = Point(s['x_wave'][i], s['y_wave'][i])

    crs_point = extract_numbers(c['crs'])
    crs_wave = extract_numbers(c['crs_wave'])

    # Ensure the input is a valid geometry object and the CRS is correctly specified
    point_df = gpd.GeoDataFrame(geometry=[wave_point], crs=f"EPSG:{crs_point[0]}")
    point_df_wgs84 = point_df.to_crs(epsg=int(crs_wave[0]))

    # Extract the coordinates for the point with index i
    point = list(point_df_wgs84.iloc[0]['geometry'].coords)[0]

    # Find the wave data file index corresponding to the closest point of the wave extraction point
    # Read the first wave data file
    first_file = glob.glob(wave_dir + r"\*.nc")[0]
    ds = xr.open_dataset(first_file, engine='netcdf4')

    geod = pyproj.Geod(ellps='WGS84')
    coord = pd.DataFrame({'x': ds.longitude[:], 'y': ds.latitude[:]})

    # Find the closest index
    lat, lon = point[1], point[0]
    distances = np.array([geod.inv(x, y, lon, lat)[2] for x, y in zip(coord['x'].values, coord['y'].values)])
    closest_index = int(np.argmin(distances))

    # Read wave data
    data_list = []

    for file in glob.glob(wave_dir + r"\*.nc"):
        ds = xr.open_dataset(file, engine='netcdf4')

        # Extract the variables (First 264 time steps are discarded - model warm-up)
        model_data_temp = pd.DataFrame({
            'Date': ds.time[264:], 
            'depth': ds.botl[264:, closest_index],
            'Hm0': ds.hs[264:, closest_index], 
            'tp': ds.tps[264:, closest_index],
            'theta': ds.theta0[264:, closest_index]
        })

        # Convert the wave direction from cartesian to nautical
        model_data_temp['theta'] = np.mod(-90 - model_data_temp['theta'], 360)

        # Append the extracted data to the list
        data_list.append(model_data_temp)

    # Concatenate the list into a single DataFrame
    combined_data = pd.concat(data_list)

    # Sort combined data by 'Date'
    combined_data.sort_values('Date', inplace=True)

    return combined_data

def extract_numbers(input_string):
    # Find all sequences of digits in the input string
    numbers = re.findall(r'\d+', input_string)
    # Convert the found sequences to integers
    return [int(num) for num in numbers]

def read_wave_data(s, c, f, i):
    #Read wave file
    wave_file = f['wave_file']
    wave_data = pd.read_csv(wave_file, sep=',')
    # Convert 'datetime' column to datetime object
    wave_data['datetime'] = pd.to_datetime(wave_data['datetime'])
    
    # Set 'datetime' as the index
    wave_data.set_index('datetime', inplace=True)
    
    # Read depth data
    depth_value = float(c['depth_wave_file'])
    wave_data['depth'] = depth_value*np.ones(len(wave_data))

    # Create a new DataFrame with the specified structure
    result_df = pd.DataFrame({
        'Date': wave_data.index,
        'depth': wave_data['depth'],
        'Hm0': wave_data['significant_wave_height_m'],
        'tp': wave_data['peak_wave_period_s'],
        'theta': wave_data['wave_direction_deg_true_north']
    }).reset_index(drop=True)

    
    return result_df

@jit(nopython=True)
def larson_direct_method(H,T,d,theta,dir_norm):
    g=9.81
    pi=np.pi
    gambr=0.78

    Hb = np.zeros(len(H))
    dirb = np.zeros(len(H))
    db = np.zeros(len(H))
    
    for i in range(len(H)):
        
        
        #Calculate angle between incoming wave and shore normal
        dir_deg=dir_norm-theta[i]

        if dir_deg > 90 or dir_deg < -90 or T[i] < 0.1:
            Hb[i] =0
            dirb[i] = 0
            db[i] = 0
        else:

            #Convert wave direction to radians
            dir = np.radians(dir_deg)

            #Determine wave properties at input location using a Pade approximation
            y=(2*pi/T[i])**2*d[i]/g
            F=y+1/(1+y*(0.66667+y*(0.3555+y*(0.16084+y*(0.06320+\
                        y*(0.02174+y*(0.00654+y*(0.00171+y*(0.00039+y*0.000111)))))))))
            L=2*pi/np.sqrt(y*F/d[i]**2)

            C=4*pi*d[i]/L

            if C < 15:
                CN=0.5*(1+C/np.sinh(C)) 
            else:
                CN=0.5

            CG=CN*L/T[i] #Group celerity

            CIN=L/T[i] #Phase celerity

            ALFA=(CIN/np.sqrt(g*H[i]))**4*gambr**2/CN

            XA=(np.cos(dir)/ALFA)**0.4

            E=(np.sin(dir))**2*XA

            DELTA=1.0+0.1649*E+0.5948*E**2-1.6787*E**3+2.8573*E**4

            X=XA*DELTA

            db[i]=X*CIN**2/g

            Hb[i]=db[i]*gambr

            if np.abs(np.sin(dir)*np.sqrt(X)) > 1:
                Hb[i] = 0
                dirb[i] = 0
            else:
                dirb[i] = np.arcsin(np.sin(dir)*np.sqrt(X))

    return Hb, dirb, db


def year_max_total_WL(df):
    # Create a DataFrame with the date and level data
    df['Year'] = pd.to_datetime(df['Date']).dt.year

    # Shift year boundary to July
    df['Year'] = df['Date'].apply(lambda x: x.year if x.month < 7 else x.year + 1)

    # Calculate the maximum water level for each year
    year_max_level = df.groupby('Year')['TWL'].max().reset_index()

    return year_max_level['TWL']

def plot_GEV_TWL(year_max_TWL, shape, loc, scale, i,s,c):
    
    #Generate a range of values for the x-axis
    x = np.linspace(min(year_max_TWL), max(year_max_TWL), num=100)
    #Calculate the probability density function (PDF) using the fitted parameters
    pdf = genextreme.pdf(x, shape, loc, scale)
    #Plot the histogram of the data and the fitted GEV PDF
    plt.figure()
    plt.hist(year_max_TWL, bins=20, density=True, alpha=0.6, label="Data Histogram")
    plt.plot(x, pdf, 'r', label="Fitted GEV PDF")
    plt.xlabel("Yearly maximum total water level [m rel MSL]")
    plt.ylabel("Probability Density")
    plt.title("Fitting GEV Distribution")
    plt.legend()

    # Save the plot to the 'out' folder
    plt.savefig(f'out/gev_distribution_TWL_{s["ID"][i]}.png')
    plt.close()

    # Calculate the return levels for different return periods
    return_periods_eq = np.linspace(0, 100, num=100)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        return_levels = genextreme.ppf(1 - 1 / np.array(return_periods_eq), shape, loc, scale)
    # Observed data
    levels = year_max_TWL

    # Step 1: Calculate the ranks
    ranks = np.argsort(np.argsort(-levels)) + 1

    # Step 2: Calculate the return periods
    N = len(levels)
    return_periods = (N + 1) / ranks


    # Plot the return levels
    plt.figure(figsize=(3.5, 3))
    #Plot the observed data
    plt.plot(return_periods, levels, 'go')
    plt.plot(return_periods_eq, return_levels, 'k-')
    plt.plot(return_periods_eq, return_levels+float(c['SLR']), 'k--')
    plt.legend(['Simulated data', rf"Fitted GEV ($\mu = {loc:.2f}$, $\sigma = {scale:.2f}$, $\xi = {shape:.2f}$)", 'Fitted GEV + SLR'], loc='upper center', bbox_to_anchor=(0.45, -0.20), ncol=1)
    plt.xlabel("Return period (years)")
    plt.ylabel("Return level (m rel MSL)")
    plt.xlim([0, 100])
    plt.title("Total water level")
    # Save the plot to the 'out' folder
    plt.savefig(f'out/return_levels_TWL{s["ID"][i]}.png', dpi=300, bbox_inches='tight')
    plt.close()


# Reverse shoaling to find the equivalent wave height at deep water
@jit(nopython=True)
def rev_shoal(T,db,Hb,dirb):
    g=9.81
    pi=np.pi

    Hprim_onshore = np.zeros(len(T))

    for i in range(len(T)):

        cgdeep=g*T[i]/(4*pi)
        cgbreak=np.sqrt(g*db[i])

        if cgbreak == 0:
            Ks = 0
        else:
            Ks=np.sqrt(cgdeep/cgbreak)

        if Ks == 0:
            Hprim = 0
        else:
            Hprim = Hb[i]/Ks

        # Correct wave height to only account for onshore component at breaking

        Hprim_onshore[i] = Hprim

#        Hprim_onshore[i] = Hprim*np.cos(dirb[i])

    return Hprim_onshore

@jit(nopython=True)
def stockdon(Hprim_onshore, T, alpha, SWL):
    g = 9.81
    pi = np.pi
    TWL = np.zeros(len(Hprim_onshore))

    for i in range(len(Hprim_onshore)):
        if T[i] == 0:
            L = 0
        else:
            L = g * T[i]**2 / (2 * pi)
        
        # Iribarren number
        if L == 0 or Hprim_onshore[i] == 0:
            ir = 0
        else:
            ir = alpha / np.sqrt(Hprim_onshore[i] / L)
        
        # Stockdon formula
        if Hprim_onshore[i] <= 0.1:
            R = 0
        elif ir < 0.3:
            R = 0.043 * np.sqrt(Hprim_onshore[i] * L)
        else:
            R = 1.1 * (0.35 * alpha * np.sqrt(Hprim_onshore[i] * L) + 0.5 * np.sqrt(Hprim_onshore[i] * L * (0.563 * alpha**2 + 0.004)))
        
        TWL[i] = R + SWL[i]
    
    return TWL


@jit(nopython=True)
def storm_erosion(zb,Hprim_onshore,T,SWL,c_impact,dt,dist,ave_height,slope,L_veg,V):
    g=9.81
    cf=0.02
    cf_veg=0.04
    pi=np.pi
    storm_erosion = np.zeros(len(Hprim_onshore))
    #Initialize xt
    distin = dist
    
    for i in range(len(Hprim_onshore)):
        #Calculate runup Larson eq
        if T[i] == 0:
            L = 0
            R = 0
        else:
            #Deep water wave length
            L=g*T[i]**2/(2*pi)
            #Convert sign wave height to Hrms
            H=Hprim_onshore[i]*0.707
            R = 0.158*np.sqrt(H*L)

        if R + SWL[i] > zb:
            #Calculate xt
            if V/L_veg < 0.2:
                xt = L_veg
            else:
                xt = dist - SWL[i]/slope

            # Adjust runup to account for friction
            if V/L_veg < 0.2:
                RPRIM = R*np.exp(-2*cf_veg*xt)+(zb-SWL[i])*(1-np.exp(-2*cf_veg*xt))
            else:
                RPRIM = R*np.exp(-2*cf*xt)+(zb-SWL[i])*(1-np.exp(-2*cf*xt))

            if RPRIM + SWL[i] > zb:
                #Calculate erosion [m3 during time step]
                storm_erosion[i] = 4*c_impact*(RPRIM+SWL[i]-zb)**2/T[i]*dt[i]
                #Update distance between the shoreline and vegetation line
                if ave_height > 0:
                    dist += storm_erosion[i]/ave_height
                else:
                    dist = L_veg
            
            else:
                storm_erosion[i] = 0
                #reset distance between shoreline and vegetation line
                dist = distin
        else:
            storm_erosion[i] = 0
            #reset distance between shoreline and vegetation line
            dist = distin
            
    return storm_erosion    

def sum_erosion_per_year(df):
    # Create a DataFrame with the date and level data
    df['Year'] = pd.to_datetime(df['Date']).dt.year

    # Shift year boundary to July
    df['Year'] = df['Date'].apply(lambda x: x.year if x.month < 7 else x.year + 1)

    # Calculate the sum of erosion for each year
    yearly_erosion = df.groupby('Year')['erosion'].sum().reset_index()

    return yearly_erosion['erosion']

def plot_GPD_erosion(yearly_erosion,excesses, shape, loc, scale, i):
        
        #Generate a range of values for the x-axis
        x = np.linspace(min(yearly_erosion), max(yearly_erosion), num=100)
        #Calculate the probability density function (PDF) using the fitted parameters
        pdf = genpareto.pdf(x, shape, loc, scale)
        #Plot the histogram of the data and the fitted GPD PDF
        plt.figure()
        plt.hist(yearly_erosion, bins=20, density=True, alpha=0.6, label="Data Histogram")
        plt.plot(x, pdf, 'r', label="Fitted GPD PDF")
        plt.xlabel("Annual Erosion Volume [m3]")
        plt.ylabel("Probability Density")
        plt.title("Fitting GPD Distribution")
        plt.legend()
    
        # Save the plot to the 'out' folder
        plt.savefig(f'out/gpd_distribution_erosion_{i}.png')
        plt.close()
    
        # Calculate the return levels for different return periods
        return_periods = [1, 2, 5, 10, 20, 50, 100]
        num_exceedances = len(excesses)
        total_data_points = len(yearly_erosion)
        exceedance_prob = num_exceedances / total_data_points
        
        # Calculate the return levels for the specified return periods
        return_levels = genpareto.ppf(1 - exceedance_prob / np.array(return_periods), shape, loc, scale)
        # Plot the return levels
        plt.figure()
        plt.plot(return_periods, return_levels, 'bo-')
        plt.xlabel("Return Period (years)")
        plt.ylabel("Return Level Erosion [m3]")
        plt.title("Return Levels for Different Return Periods")
        # Save the plot to the 'out' folder
        plt.savefig(f'out/return_levels_erosion_{i}.png')
        plt.close()