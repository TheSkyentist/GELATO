"""Utility Functions for GELATO"""

# Import packages
import sys
import argparse
import numpy as np
from os import path
from datetime import datetime
from astropy.io import fits
from astropy.table import Table

def header():
    print("Welcome to GELATO")
    print("Galaxy/AGN Emission Line Analysis TOol")
    print("Developed by R. E. Hviding")
    now = datetime.now()
    print("Started making gelato at",now)
    return now

def footer(then):
    now = datetime.now()
    print("Finished making GELATO at",datetime.now())
    print("Elapsed time:",now - then)

def parseArgs():
    
    # Parse Arguements
    parser = argparse.ArgumentParser()
    parser.add_argument('Parameters', type=str, help='Path to parameters file.')
    parser.add_argument('-s','--single', action='store_true', help='Specify single object mode.')
    if ((len(sys.argv) > 2) and (sys.argv[2][0] == '-')):
        parser.add_argument('Spectrum', type=str, help='Path to spectrum.')
        parser.add_argument('Redshift', type=float, help='Redshift of object.')
    else:
        parser.add_argument('ObjectList', type=str, help='Path to object list with paths to spectra and their redshifts.')
    return parser.parse_args()

def loadObjects(tpath):
    if tpath.endswith('.csv'):
        objects = np.atleast_1d(np.genfromtxt(tpath,delimiter=',',dtype=['U100',np.float_],names=['Path','z']))
    elif tpath.endswith('.fits'):
        objects = Table.read(tpath)
        objects.convert_bytestring_to_unicode()
        objects = np.atleast_1d(objects)
    else:
        print('Object list not .csv or .fits.')
        exit()
    return objects

# Get filename
def fileName(fpath):

    # Supported file types
    ftypes = ['.fit','.fits','.fit.gz','fits.gz']

    # If a supported filetype, return w/o the filetype
    for ft in ftypes:
        if (fpath[-len(ft):] == ft): return fpath[:-len(ft)]

    # Otherwise just return it
    return fpath

# Loading Bar (Deprecated)
def loadingBar(i,N,L=40):

    # If at 100%, print the final part of the bar
    if i == N:
        print('Progress: |'+L*'#'+'| 100%')
        return

    # Loading bar params
    p = int(100*i/N) # Percentage
    l = int(L*i/N) # Length of bar
    if p == 0:
        print('Progress: |'+L*'-'+'|   0%',end='\r')
    elif p < 10:
        print('Progress: |'+l*'#'+(L-l)*'-'+'|   '+str(p)+'%',end='\r')
    else:
        print('Progress: |'+l*'#'+(L-l)*'-'+'|  '+str(p)+'%',end='\r')