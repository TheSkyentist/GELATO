#! /usr/bin/env python

""" Wrapper for mulitple GELATO runs """

# Packages
import os
import copy
import argparse
import numpy as np

# GELATO supporting files
import GELATO
import ConstructParams as CP

# Main Function
if __name__ == "__main__":

    ## Parse Arguements to find Parameter File ##
    parser = argparse.ArgumentParser()
    parser.add_argument('Parameters', type=str, help='Path to parameters file')
    parser.add_argument('ObjectList', type=str, help='Path to object list with paths to spectra and their redshifts.')
    args = parser.parse_args()
    p = CP.construct(args.Parameters)
    ## Parse Arguements to find Parameter File ##

    ## Create Directory for Output
    if not os.path.exists(p["OutFolder"]):
        os.mkdir(p["OutFolder"])

    if p['Verbose']:
        GELATO.header()

    ## Assemble Objects
    objects = np.genfromtxt(args.ObjectList,delimiter=',',dtype=['U100',np.float_],names=['File','z'])
    ## Assemble Objects

    ## Run GELATO ##
    if p['NProcess'] > 1: # Mutlithread
        import multiprocessing as mp
        pool = mp.Pool(processes=p['NProcess'])
        inputs = [(copy.deepcopy(p),o['File'],o['z']) for o in objects]
        pool.starmap(GELATO.GELATO, inputs)
        pool.close()
        pool.join()
    else: # Single Thread
        for o in objects: GELATO.GELATO(copy.deepcopy(p),o['File'],o['z'])
    ## Run GELATO ##

    ## Concatenate Results ##
    if p['Concatenate']:
        import ConcatResults as CR
        CR.concatfromresults(p,objects)
    ## Concatenate Results ##

    if p['Verbose']:
        GELATO.footer()