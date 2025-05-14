# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:28:58 2023

@author: Caroline Hallin

This script reads a configuration file and updates the dictionaries c, f, and p based on the key-value pairs found in the file.
The configuration file is expected to be in a specific format, where each line contains a key-value pair separated by an equal sign (=).
"""

import os
from dictionaries import *


# Specify the relative directory path where the input file is located
dirname = os.path.dirname(__file__)
input_path = os.path.join(dirname, 'in')


# Specify the name of the input file (LATER: ADD THE POSSIBILITY TO DEFINE THE FILE NAME IN THE COMMAND LINE)
input_file_name = 'input.txt'

# Check if the input file exists in the specified directory
if not os.path.exists(input_path):
    print(f"The directory '{input_path}' does not exist.")
    exit()

if not os.path.isfile(os.path.join(input_path, input_file_name)):
    print(f"The file '{input_file_name}' does not exist in the specified directory.")
    exit()

# Read the input file
open(os.path.join(input_path, input_file_name), 'r')
with open(os.path.join(input_path, input_file_name), 'r') as file:
    for line in file:
        line = line.strip()
        # Remove empty lines and lines starting with % or #
        if not line or line.startswith("%") or line.startswith("#"):
            continue
        # Remove whitespaces and text after hashtags
        cleaned_line = line.split('#')[0].strip()
        # Split the line into key and new value
        key, value = cleaned_line.strip().split('=')
        key = key.strip()
        value = value.strip()
        
        # Update the value in the existing c dictionary if the key exists
        if key in c:
            c[key] = value

        # Update the value in the existing f dictionary if the key exists
        if key in f:
            f[key] = os.path.join(input_path, value)

        # Update the value in the existing p dictionary if the key exists
        if key in p:
            p[key] = value.lower() == 'true'

