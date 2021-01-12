#! /usr/bin/env python

""" Wrapper for single GELATO run """

# Packages
import os
import sys
import argparse

# GELATO supporting files
import GELATO
import ConstructParams as CP

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
        os.mkdir(p["OutFolder"])

    if p['Verbose']:
        GELATO.header()

    ## Run GELATO ##
    GELATO.GELATO(p, args.Spectrum, args.Redshift)

    if p['Verbose']:
        GELATO.footer()