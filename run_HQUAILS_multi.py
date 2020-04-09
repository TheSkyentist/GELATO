#! /usr/bin/env python

""" Wrapper for mulitple HQUAILS runs """

# Packages
import copy
import argparse
import numpy as np

# HQUAILS supporting files
import HQUAILS
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

    if p['Verbose']:
        HQUAILS.header()


    ## Assemble Objects
    objects = np.genfromtxt(args.ObjectList,delimiter=',',dtype=['U100',np.float_],names=['File','z'])
    ## Assemble Objects

    ## Run HQUAILS ##
    if p['NProcess'] > 1: # Mutlithread
        import multiprocessing as mp
        pool = mp.Pool(processes=p['NProcess'])
        inputs = [(copy.deepcopy(p),o['File'],o['z']) for o in objects]
        pool.starmap(HQUAILS.HQUAILS, inputs)
        pool.close()
        pool.join()
    else: # Single Thread
        for o in objects: HQUAILS.HQUAILS(copy.deepcopy(p),o['File'],o['z'])
    ## Run HQUAILS ##

    ## Concatenate Results ##
    if p['Concatenate']:
        import ConcatResults as CR
        CR.concatfromresults(p,objects)
    ## Concatenate Results ##

    if p['Verbose']:
        HQUAILS.footer()