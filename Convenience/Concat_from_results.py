#! /usr/bin/env python

""" Convinience Function to Concatenate Results """

# Packages
import argparse
import numpy as np
from astropy.table import Table

# GELATO
import gelato.Concatenate as C
import gelato.ConstructParams as CP

## Parse Arguements ##
parser = argparse.ArgumentParser()
parser.add_argument('Parameters', type=str, help='Path to parameters file')
parser.add_argument('ObjectList', type=str, help='Path to object list with paths to spectra and their redshifts.')
args = parser.parse_args()
p = CP.construct(args.Parameters)
## Parse Arguements

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

## Concatenate Results
C.concatfromresults(p,objects)