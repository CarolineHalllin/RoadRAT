# -*- coding: utf-8 -*-
"""
RoadRAT - a framework to calculate the probability of inundation, runup, and storm erosion impacting coastal roads 
in present and future climate scenarios.

This is the main program for RoadRAT. It imports all necessary functions and modules, sets up the computation points, 
and runs the calculations for inundation, runup, and erosion.
It also handles the output of results and plotting of maps.

Caroline Hallin 2023-08-29
"""

# Import functions from other scripts
from dictionaries import *
from input import *
from longterm import *
from runup_erosion import *
from SWL import *
from transects import *
from output import *
import logging
import pprint
import os

# Set up logging
if not os.path.exists('logs'):
    os.makedirs('logs')

logging.basicConfig(
    filename='logs/road_rat.log',
    filemode='w',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Log to console
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_formatter = logging.Formatter('%(levelname)s - %(message)s')
console_handler.setFormatter(console_formatter)
logging.getLogger().addHandler(console_handler)


def main():

    # Check if output directory exists, if not create it
    if not os.path.exists('out'):
        os.makedirs('out')

    #Log dictionaries
    logging.info("Updated dictionaries:\nc:\n%s\nf:\n%s\np:\n%s",
             pprint.pformat(c),
             pprint.pformat(f),
             pprint.pformat(p))
     
    # Extract computation points and populate with input data
    computation_points(s,c,f)

    # Print output-file
    output(s,c['output_file'])

    # Calculate probability of flooding due to SWL exceeding road level
    if p['process_inundation']:
        GEV_SWL(s,c,f,p)

    # Print output-file
    output(s,c['output_file'])
        
    # Derive morphological parameters for runup and erosion calculations
    if p['process_runup'] or p['process_erosion']:

        morph_parameters(s,c,f,p)

        # Print output-file
        output(s,c['output_file'])

        # Calculate probability of runup or erosion reaching the road
        prob_runup_erosion(s,c,f,p)

        # Print output-file
        output(s,c['output_file'])
        

    if p['process_plotting']:
        # Plot map with output
        plot_results(s,c,f,p)

    # Print netcdf output file
    output_netcdf(s,c['output_file'])




main()    

