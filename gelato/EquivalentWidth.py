#! /usr/bin/env python

""" Generate Rest Equivalent Widths """

# Import packages
import numpy as np
from astropy.io import fits
from astropy.table import Table,hstack

# gelato supporting files
import gelato.CustomModels as CM
import gelato.SpectrumClass as SC

# Calculate Equivalent Width
def EquivalentWidth(spectrum,model,param_names,parameters):

    # Continuum
    if 'PL_Continuum_Coefficient' in param_names:
        continuum = model[0:2]
    else: continuum = model[0]

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
                if params["Verbose"]:
                    print("Equivalent Widths Already Generated:",path.split('/')[-1])
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
            namesplit[-1] = 'REW'
            EWnames.append('-'.join(namesplit))
            
            # Get line center
            center = float(namesplit[-2])

            # Get line flux/redshift
            lineflux = parameters[names[i+1]]
            oneplusz = 1+parameters[names[i]]

            # Initialize continuum fluxes
            contflux = np.ones(len(parameters))

            # Iterate over params
            for j,params in enumerate(parameters):

                # Load in continuum
                continuum.parameters = params[:ind]
                contflux[j] = continuum(center*oneplusz[j])

            # Get EW
            EWs.append(np.abs(lineflux/(contflux*oneplusz)))

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
    elif single: # One EW
        EWfromresults(p, args.Spectrum, args.Redshift)
    elif multi: # Many EW
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