#! /usr/bin/env python

""" Convinience Function to Generate Rest Equivalent Widths """

# Packages
import sys
import copy
import argparse
import numpy as np
from astropy.table import Table

# GELATO
import gelato.ConstructParams as CP
import gelato.EquivalentWidth as EW

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
elif single: # One EW
    EW.EWfromresults(p, args.Spectrum, args.Redshift)
elif multi: # Many EW
    ## Assemble Objects
    if args.ObjectList.endswith('.csv'):
        objects = np.atleast_1d(np.genfromtxt(args.ObjectList,delimiter=',',dtype=['U100',np.float_],names=['Path','z']))
    elif args.ObjectList.endswith('.fits'):
        objects = Table.read(args.ObjectList)
        objects.convert_bytestring_to_unicode()
        objects = np.atleast_1d(objects)
    else:
        print('Object list not .csv or .fits.')
    ## Assemble Objects
    if p['NProcess'] > 1: # Mutlithread
        import multiprocessing as mp
        pool = mp.Pool(processes=p['NProcess'])
        inputs = [(copy.deepcopy(p),o['Path'],o['z']) for o in objects]
        pool.starmap(EW.Wfromresults, inputs)
        pool.close()
        pool.join()
    else: # Single Thread
        for o in objects: EW.EWfromresults(copy.deepcopy(p),o['Path'],o['z'])