#! /usr/bin/env python

""" Convinience Function for Plotting """

# Packages
import sys
import copy
import argparse
import numpy as np

# GELATO
import gelato.Plotting as P
import gelato.ConstructParams as CP

## Parse Arguements to find Parameter File ##
parser = argparse.ArgumentParser()
parser.add_argument('Parameters', type=str, help='Path to parameters file')
parser.add_argument('--ObjectList', type=str, help='Path to object list with paths to spectra and their redshifts.')
parser.add_argument('--Spectrum', type=str, help='Path to spectrum.')
parser.add_argument('--Redshift', type=float, help='Redshift of object')
args = parser.parse_args()
p = CP.construct(args.Parameters)
## Parse Arguements to find Parameter File ##

# Check if we are doing single or multi
single = args.Spectrum != None and args.Redshift != None
multi = args.ObjectList != None

if single == multi:
    print('Specify either Object List XOR Spectrum and Redshift.')
    print('Both or neither were entered.')
elif single: # One Plot
    P.plotfromresults(p, args.Spectrum, args.Redshift)
elif multi: # Many plots
    # Load Objects
    objects = np.genfromtxt(args.ObjectList,delimiter=',',dtype='U100,f8',names=['File','z'])
    if p['NProcess'] > 1: # Mutlithread
        import multiprocessing as mp
        pool = mp.Pool(processes=p['NProcess'])
        inputs = [(copy.deepcopy(p),o['File'],o['z']) for o in objects]
        pool.starmap(P.plotfromresults, inputs)
        pool.close()
        pool.join()
    else: # Single Thread
        for o in objects: P.plotfromresults(copy.deepcopy(p),o['File'],o['z'])