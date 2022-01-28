#! /usr/bin/env python

""" Convinience Function to Generate Rest Equivalent Widths """

# Packages
import copy

# GELATO
import gelato.Utility as U
import gelato.ConstructParams as CP
import gelato.EquivalentWidth as EW

# Main Function
if __name__ == "__main__":

    # Parse Arguement
    args = U.parseArgs()

    # Parameters
    p = CP.construct(args.Parameters)

    # Single Mode
    if args.single:

        EW.EWfromresults(p, args.Spectrum, args.Redshift)

    # Multi Mode
    else:

        ## Assemble Objects
        objects = U.loadObjects(args.ObjectList)

        # EWs
        if p['NProcess'] > 1: # Mutlithread
            import multiprocessing as mp
            pool = mp.Pool(processes=p['NProcess'])
            inputs = [(copy.deepcopy(p),o['Path'],o['z']) for o in objects]
            pool.starmap(EW.EWfromresults, inputs)
            pool.close()
            pool.join()
        else: # Single Thread
            for o in objects: EW.EWfromresults(copy.deepcopy(p),o['Path'],o['z'])