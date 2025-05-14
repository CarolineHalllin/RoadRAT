# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 13:57:09 2023

@author: Caroline Hallin

Script to write output files and plot results.
"""

import os
import numpy as np
import rasterio
from dictionaries import *
import matplotlib.pyplot as plt
import shapefile
import xarray as xr
import logging
from shapely.geometry import LineString
from matplotlib.colors import ListedColormap, BoundaryNorm


def output(s, filename):
    filename = os.path.join('out', filename)

    # Define the header based on the keys in the dictionary `s`
    header = list(s.keys())

    # Determine the maximum width for each column
    col_widths = {
        key: max(
            len(str(key)),
            max((len(f"{item:.2f}") if isinstance(item, float) else len(str(item))) for item in s[key]) if len(s[key]) > 0 else len(str(key))
        )
        for key in header
    }

    # Open file for writing
    with open(filename, 'w') as file:
        logging.info(f'Writing text output to {filename}')
        # Write the header with adjusted widths
        file.write('\t'.join(f"{key:<{col_widths[key]}}" for key in header) + '\n')
        
        # Determine the number of rows by finding the maximum length of the nested dictionaries
        num_rows = max(len(v) for v in s.values())

        # Write each row of data
        for i in range(num_rows):
            row = []
            for key in header:
                if i < len(s[key]):
                    value = s[key][i]
                    if isinstance(value, float):
                        formatted_value = f"{value:.2f}"
                    else:
                        formatted_value = str(value)
                else:
                    formatted_value = ''
                row.append(f"{formatted_value:<{col_widths[key]}}")
            file.write('\t'.join(row) + '\n')

def output_netcdf(s, filename):
    logging.info(f'Writing NetCDF output to {filename}')
    filename = os.path.join('out', filename)
    #replace the file extension with .nc
    filename = os.path.splitext(filename)[0] + '.nc'

    # Create a new NetCDF file
    with xr.Dataset() as ds:
        # Add variables to the dataset
        for key in s.keys():
            ds[key] = xr.DataArray(s[key], dims=['points'], coords={'points': range(len(s[key]))})
            

        # Save the dataset to a NetCDF file
        ds.to_netcdf(filename, mode='w')


def plot_results(s, c, f, p):
    # Plot the results of the calculations
    params_to_plot = []
    if 'p_swl' in s and len(s['p_swl']) > 0:
        s['p_swl'] = np.array(s['p_swl'])
        s['p_swl'][np.isinf(s['p_swl'])] = 1000
        params_to_plot.append('p_swl')
    if 'p_r' in s and len(s['p_r']) > 0:
        s['p_r'] = np.array(s['p_r'])
        s['p_r'][np.isinf(s['p_r'])] = 1000
        params_to_plot.append('p_r')
    if 'p_er' in s and len(s['p_er']) > 0:
        s['p_er'] = np.array(s['p_er'])
        s['p_er'][np.isinf(s['p_er'])] = 1000
        params_to_plot.append('p_er')
    if 'p_swl_slr' in s and len(s['p_swl_slr']) > 0:
        s['p_swl_slr'] = np.array(s['p_swl_slr'])
        s['p_swl_slr'][np.isinf(s['p_swl_slr'])] = 1000
        params_to_plot.append('p_swl_slr')
    if 'p_r_slr' in s and len(s['p_r_slr']) > 0:
        s['p_r_slr'] = np.array(s['p_r_slr'])
        s['p_r_slr'][np.isinf(s['p_r_slr'])] = 1000
        params_to_plot.append('p_r_slr')
    if 'p_er_slr' in s and len(s['p_er_slr']) > 0:
        s['p_er_slr'] = np.array(s['p_er_slr'])
        s['p_er_slr'][np.isinf(s['p_er_slr'])] = 1000
        params_to_plot.append('p_er_slr')

    categories = [0, 2, 10, 20, 40, 60, 80, 100]
    colors = ['darkred', 'red', 'orange', 'yellow', 'green', 'cyan', 'blue']
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(categories, cmap.N)

    sf = shapefile.Reader(f['road_file'])
    shape = sf.shape(0)
    points = shape.points
    line = LineString(points)

    with rasterio.open(f['DEM_file']) as dem:
        dem_data = dem.read(1)
        dem_extent = (dem.bounds.left, dem.bounds.right, dem.bounds.bottom, dem.bounds.top)
        z_min = 0
        z_max = max(s['z']) + 1

        for param in params_to_plot:
            fig, ax = plt.subplots(figsize=(10, 5))
            fig.suptitle('RoadRAT results')

            ax.imshow(dem_data, extent=dem_extent, cmap='terrain', origin='upper', vmin=z_min, vmax=z_max)
            ax.plot(*line.xy, color='black', linewidth=2)
            ax.set_aspect('equal', adjustable='box')

            exceed_mask = np.array(s[param]) > categories[-1]
            within_mask = ~exceed_mask

            s['x'] = np.array(s['x'])
            s['y'] = np.array(s['y'])

            if np.any(within_mask):
                scatter = ax.scatter(s['x'][within_mask], s['y'][within_mask], c=np.array(s[param])[within_mask], cmap=cmap, norm=norm, edgecolors='black', linewidths=1, s=100)

            if np.any(exceed_mask):
                ax.scatter(s['x'][exceed_mask], s['y'][exceed_mask], facecolors='none', edgecolors='black', linewidths=1, s=100)

            if param == 'p_swl':
                ax.set_title(f'Return period inundation SWL ({c["start_yr"]})')
            elif param == 'p_r':
                ax.set_title(f'Return period exposure to runup ({c["start_yr"]})')
            elif param == 'p_er':
                ax.set_title(f'Return period damage due to erosion ({c["start_yr"]})')
            elif param == 'p_swl_slr':
                ax.set_title(f'Return period inundation SWL with SLR ({c["end_yr"]})')
            elif param == 'p_r_slr':
                ax.set_title(f'Return period exposure to runup with SLR ({c["end_yr"]})')
            elif param == 'p_er_slr':
                ax.set_title(f'Return period damage due to erosion with SLR ({c["end_yr"]})')

            # Zoom in on the extent of the road segment
            minx, miny, maxx, maxy = line.bounds
            buffer = 2 * float(c['spacing'])
            plt.xlim(minx - buffer, maxx + buffer)
            plt.ylim(miny - buffer, maxy + buffer)

            if 'scatter' in locals():
                cbar = fig.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.1, aspect=50)
                cbar.set_label('Return Period (years)')
                # Calculate the midpoints between the category boundaries
                tick_positions = [(categories[i] + categories[i+1]) / 2 for i in range(len(categories) - 1)]
    
                cbar.set_ticks(tick_positions)
                cbar.set_ticklabels(['0-2', '2-10', '10-20', '20-40', '40-60', '60-80', '80-100'])

            plt.savefig(f'out/RoadRAT_results_{param}.png')
            plt.close(fig)