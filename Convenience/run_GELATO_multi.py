#! /usr/bin/env python

""" Wrapper for mulitple gelato runs """

# Packages
import os
import copy
import argparse
import numpy as np

# gelato supporting files
import gelato
import gelato.ConstructParams as CP

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
        gelato.header()

    ## Assemble Objects
    objects = np.genfromtxt(args.ObjectList,delimiter=',',dtype=['U100',np.float_],names=['File','z'])
    ## Assemble Objects

    ## Run gelato ##
    if p['NProcess'] > 1: # Mutlithread
        import multiprocessing as mp
        pool = mp.Pool(processes=p['NProcess'])
        inputs = [(copy.deepcopy(args.Parameters),o['File'],o['z']) for o in objects]
        pool.starmap(gelato.gelato, inputs)
        pool.close()
        pool.join()
    else: # Single Thread
        for o in objects: gelato.gelato(copy.deepcopy(args.Parameters),o['File'],o['z'])
    ## Run gelato ##

    ## Concatenate Results ##
    if p['Concatenate']:
        import gelato.ConcatResults as CR
        CR.concatfromresults(p,objects)
    ## Concatenate Results ##

    if p['Verbose']:
        gelato.footer()