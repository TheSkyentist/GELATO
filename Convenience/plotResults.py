#! /usr/bin/env python

""" Convinience Function for Plotting """

# Packages
import copy

# GELATO
import gelato.Utility as U
import gelato.Plotting as P
import gelato.ConstructParams as CP

# Main Function
if __name__ == "__main__":

    # Parse Arguement
    args = U.parseArgs()

    # Parameters
    p = CP.construct(args.Parameters)

    # Single Mode
    if args.single:

        P.plotfromresults(p, args.Spectrum, args.Redshift)

    # Multi Mode
    else:

        # Assemble Objects
        objects = U.loadObjects(args.ObjectList)

        # Plot
        if p['NProcess'] > 1: # Mutlithread
            import multiprocessing as mp
            pool = mp.Pool(processes=p['NProcess'])
            inputs = [(copy.deepcopy(p),o['Path'],o['z']) for o in objects]
            pool.starmap(P.plotfromresults, inputs)
            pool.close()
            pool.join()
        else: # Single Thread
            for o in objects: P.plotfromresults(copy.deepcopy(p),o['Path'],o['z'])