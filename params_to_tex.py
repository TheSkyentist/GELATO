#! /usr/bin/env python

""" Turn parameter file into LaTeX table """

# Packages
import os
import copy
import argparse
import numpy as np

# gelato supporting files
import gelato
import gelato.ConstructParams as CP

# Main Function
if __name__ == "__main__":

    ## Parse Arguements to find Parameter File ##
    parser = argparse.ArgumentParser()
    parser.add_argument('Parameters', type=str, help='Path to parameters file')
    args = parser.parse_args()
    p = CP.construct(args.Parameters)
    ## Parse Arguements to find Parameter File ##

    ## Create Directory for Output
    if '.json' in args.Parameters:
        outname = args.Parameters.replace('.json','.tex')
    else: 
        outname = args.Parameters + '.tex'

    