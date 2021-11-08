""" Additional components """

import numpy as np
import gelato.CustomModels as CM
from gelato.Constants import C

def ComponentName(index):

    '''
    For each bit position, from right to left (increasing order of magnitude), what is the name given to each additional component such that they can be identified in the list of model name parameters.
    So the rightmost bit corresponds to index 0.
    Can be as many as desired as long as paired with a component below. 
    '''

    if index == 0:
        return 'Broad'

    if index == 1:
        return 'Outflow'

    if index == 2:
        return 'Absorption'

def AddComponent(flag, line, spectrum, prefix):

    '''
    For each bit position, what is the model that we implement,
    and how do we set the initial parameters.
    Must return astropy model with edits
    '''
    
    # Find index
    flag = bin(flag)[3:]
    index = len(flag) - flag.index('1') - 1

    ''' Broad '''
    if index == 0: 
    
        # Broad line model
        model = CM.SpectralFeature(center = line,spec = spectrum, prefix=prefix)
        x0 = model.starting()
        
        # Use wider default dispersion
        x0[2] = 1200 # Dispersion

        # Reassign Dispersion bounds
        oldMaxDispersion = model.Dispersion_bounds[1]
        model.Dispersion_bounds = (1000,10000)

        # Reassign Flux Bounds
        Fbound = 1.5*x0[1]*model.Dispersion_bounds[1]/oldMaxDispersion
        model.Flux_bounds = (-Fbound,Fbound)
        
    ''' Outflow '''
    if index == 1: 
    
        # Outflow model
        model = CM.SpectralFeature(center = line,spec = spectrum, prefix=prefix)
        x0 = model.starting()

        # Use wider default dispersion
        x0[2] = 500 # Dispersion

        # Reassign Dispersion bounds
        model.Dispersion_bounds = (500,1000)

        # Reassign Redshift bounds and starting guess
        model.Redshift_bounds = (C*(spectrum.z-0.001),C*spectrum.z)
        x0[0] = C*(spectrum.z-0.0005)
        

    return model,x0