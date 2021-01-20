""" Additional components """

import numpy as np
import gelato.CustomModels as CM

SQRT_2_PI = np.sqrt(2*np.pi)

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

def AddComponent(flag, line, spectrum):

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
        # Use wider default dispersion
        model = CM.SpectralFeature(center = line,spectrum = spectrum, Dispersion = 1000)

        # Reassignflux bounds
        oldMaxDispersion = model.Dispersion.bounds[1]
        model.Dispersion.bounds = (750,10000)
        model.Flux.bounds = (0,1.5*model.Flux*model.Dispersion.bounds[1]/oldMaxDispersion)
        
        return model

    ''' Outflow '''
    if index == 1: 
    
        # Broad line model
        # Use wider default dispersion
        model = CM.SpectralFeature(center = line,spectrum = spectrum, Dispersion = 500)
        model.Dispersion.bounds = (500,1000)

        return model

    ''' Absorption '''
    if index == 2: 
    
        # Absorption line model
        # Use wider default dispersion        
        model = CM.SpectralFeature(center = line,spectrum = spectrum, Dispersion = 600)

        # Reassign Flux
        oldMaxDispersion = model.Dispersion.bounds[1]
        oldFlux = model.Flux.value
        model.Dispersion.bounds = (350,3000)
        model.Flux.bounds = (-model.Flux.bounds[1]*1000/oldMaxDispersion,0) # Must be non-positive

        model.Redshift = spectrum.z
        model.Flux = -oldFlux*model.Dispersion/oldMaxDispersion/2
        return model