#! /usr/bin/env python

""" Convinience Function to Concatenate Results """

# Packages
import argparse
import numpy as np
from astropy.table import Table

# GELATO
import gelato.Utility as U
import gelato.Concatenate as C
import gelato.ConstructParams as CP

# Main Function
if __name__ == "__main__":

    # Parse Arguement
    args = U.parseArgs()

    # Parameters
    p = CP.construct(args.Parameters)

    # Assemble Objects
    if args.single: # Single Mode
        objects = Table([[args.Spectrum],[args.Redshift]],names=('Path','z'))
    else: # Multi Mode
        objects = U.loadObjects(args.ObjectList)

    ##Concatenate Results
    C.concatfromresults(p,objects)