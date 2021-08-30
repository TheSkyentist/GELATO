#! /usr/bin/env python

""" Convinience Function to Concatenate Results """

# Packages
import argparse
import numpy as np
from astropy.table import Table

# GELATO
import gelato
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
objects = gelato.loadObjects(args.ObjectList)
## Assemble Objects

## Concatenate Results
C.concatfromresults(p,objects)