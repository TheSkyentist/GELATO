#! /usr/bin/env python

""" Wrapper for single gelato run """

# Packages
import os
import sys
import argparse
from pathlib import Path


# gelato supporting files
import gelato
import gelato.ConstructParams as CP

# Main Function
if __name__ == "__main__":

    ## Parse Arguements to find Parameter File ##
    parser = argparse.ArgumentParser()
    parser.add_argument('Parameters', type=str, help='Path to parameters file')
    parser.add_argument('Spectrum', type=str, help='Path to spectrum.')
    parser.add_argument('Redshift', type=float, help='Redshift of object')
    args = parser.parse_args()
    p = CP.construct(args.Parameters)
    ## Parse Arguements to find Parameter File ##

    ## Create Directory for Output
    if not os.path.exists(p["OutFolder"]):
        Path(p["OutFolder"]).mkdir(parents=True)

    if p['Verbose']:
        gelato.header()

    ## Run gelato ##
    gelato.gelato(args.Parameters, args.Spectrum, args.Redshift)

    if p['Verbose']:
        gelato.footer()