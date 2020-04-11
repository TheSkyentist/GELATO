#! /usr/bin/env python

""" Generate Equivalent Widths """

# Import packages
import numpy as np
import CustomModels as CM
import SpectrumClass as SC
from astropy.io import fits
from astropy.table import Table,hstack

# Plot from results
def EWfromresults(params,path,z):

    if params["Verbose"]:
        print("Generating Equivalent Widths:",path.split('/')[-1])

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(path,z,params)

    if spectrum.regions != []:

        # Load name and parameters
        fname = params['OutFolder']+path.split('/')[-1].replace('.fits','')+'-results.fits'
        parameters = fits.getdata(fname)
        names = parameters.columns.names
        
        # Dont add if already has EWs
        for n in names:
            if 'EW' in n:
                return

        # Index where emission lines begin
        ind = (params['ContinuumDeg']+1)*len(spectrum.regions)

        ## Create continuum ##
        continuum = []
        # Add continuum
        for region in spectrum.regions:
            continuum.append(CM.Continuum(params['ContinuumDeg'],region))
        continuum = np.sum(continuum)

        EWs = []
        EWnames = []
        # Iterate over emission lines
        for i in range(ind,len(names)-1,3):

            # Split name
            namesplit = names[i].split('-')
            namesplit[-1] = 'EW'
            EWnames.append('-'.join(namesplit))
            
            # Get line center
            center = float(namesplit[-2])

            # Initialize continuum fluxes
            contflux = np.ones(len(parameters))

            # Iterate over params
            for j,params in enumerate(parameters):

                # Load in continuum
                continuum.parameters = params[:ind]
                contflux[j] = continuum(center*(1+params[names[i]]))

            # Get EW
            EWs.append(np.abs(parameters[names[i+1]]/contflux))
        
        # Save
        hstack([Table(parameters),Table(data=EWs,names=EWnames)]).write(fname,overwrite=True)
    
# EW from results
# Plot from results
if __name__ == "__main__":

    # Import if we need them
    import sys
    import copy
    import argparse
    import ConstructParams as CP

    ## Parse Arguements to find Parameter File ##
    parser = argparse.ArgumentParser()
    parser.add_argument('Parameters', type=str, help='Path to parameters file')
    parser.add_argument('--ObjectList', type=str, help='Path to object list with paths to spectra and their redshifts.')
    parser.add_argument('--Spectrum', type=str, help='Path to spectrum.')
    parser.add_argument('--Redshift', type=float, help='Redshift of object')
    args = parser.parse_args()
    p = CP.construct(args.Parameters)
    ## Parse Arguements to find Parameter File ##

    # Check if we are doing single or multi
    single = args.Spectrum != None and args.Redshift != None
    multi = args.ObjectList != None

    if single == multi:
        print('Specify either Object List XOR Spectrum and Redshift.')
        print('Both or neither were entered.')
    elif single: # One Plot
        EWfromresults(p, args.Spectrum, args.Redshift)
    elif multi: # Many plots
        # Load Obkects
        objects = np.genfromtxt(args.ObjectList,delimiter=',',dtype='U100,f8',names=['File','z'])
        if p['NProcess'] > 1: # Mutlithread
            import multiprocessing as mp
            pool = mp.Pool(processes=p['NProcess'])
            inputs = [(copy.deepcopy(p),o['File'],o['z']) for o in objects]
            pool.starmap(EWfromresults, inputs)
            pool.close()
            pool.join()
        else: # Single Thread
            for o in objects: EWfromresults(copy.deepcopy(p),o['File'],o['z'])