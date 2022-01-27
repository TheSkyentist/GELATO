#! /usr/bin/env python

""" Wrapper for mulitple gelato runs """

# Packages
import os
import copy
from pathlib import Path

# GELATO
import gelato
import gelato.Utility as U
import gelato.ConstructParams as CP

# Main Function
if __name__ == "__main__":

    # Parse Arguement
    args = U.parseArgs()

    # Parameters
    p = CP.construct(args.Parameters)

    ## Create Directory for Output
    if not os.path.exists(p["OutFolder"]):
        Path(p["OutFolder"]).mkdir(parents=True)

    if p['Verbose']:
        now = U.header()

    # Single Mode
    if args.single:

        gelato.gelato(args.Parameters, args.Spectrum, args.Redshift)

    # Multi Mode
    else:

        # Assemble Objects
        objects = U.loadObjects(args.ObjectList)

        ## Run gelato ##
        if p['NProcess'] > 1: # Mutlithread
            import multiprocessing as mp
            pool = mp.Pool(processes=p['NProcess'])
            inputs = [(copy.deepcopy(args.Parameters),o['Path'],o['z']) for o in objects]
            pool.starmap(gelato.gelato, inputs)
            pool.close()
            pool.join()
        else: # Single Thread
            for o in objects:
                gelato.gelato(copy.deepcopy(args.Parameters),o['Path'],o['z'])
        ## Run gelato ##

        # Concatenate Results
        if p['Concatenate']:
            import gelato.Concatenate as C
            C.concatfromresults(p,objects)
    
    if p['Verbose']:
        U.footer(now)