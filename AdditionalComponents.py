""" Additional components """

import numpy as np
import CustomModels as CM

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
        model = CM.SpectralFeature(center = line,spectrum = spectrum, Dispersion = 500)
        model.Dispersion.bounds = (1000,10000)

        return model

    ''' Absorption '''
    if index == 1: 
    
        # Absorption line model
        # Use wider default dispersion        
        model = CM.SpectralFeature(center = line,spectrum = spectrum, Dispersion = 600)
        originalsigma = model.Dispersion.bounds[1]
        originalflux = model.Flux.value
        model.Dispersion.bounds = (350,1000)
        model.Flux.bounds = (-model.Flux.bounds[1]*1000/originalsigma,0) # Must be non-positive

        model.Redshift = spectrum.z - 0.005
        model.Flux = -originalflux*model.Dispersion/originalsigma/2
        return model