""" Fit Model for Spectrum """

# Packages
import copy
import numpy as np
from astropy.modeling import fitting
fit = fitting.LevMarLSQFitter()

# HQUAILS supporting files
import BuildModel as BD
import FischerTest as FT
import AdditionalComponents as AC

import matplotlib.pyplot as plt

# Construct Full Model with F-tests for additional parameters
def FitComponents(spectrum,base_model,base_param_names):

    # Fit first model
    base_model = FitModel(spectrum,base_model)

    # Find number of flags
    flags = 0
    for group in spectrum.p['EmissionGroups']:
        for species in group['Species']:
            flagbits = bin(species['Flag'])[2:]
            flags += len(flagbits.replace('0',''))

    # F-tests
    accepted = []
    accepted_models = []
    accepted_names = []
    for i in range(flags):
        
        # Add new component
        EmissionGroups = AddComplexity(spectrum.p['EmissionGroups'],i)
        model,param_names = BD.BuildModel(spectrum,EmissionGroups)

        # Source starting parameters
        model = SourceParams(model,param_names,base_model,base_param_names)

        # Fit model
        model = FitModel(spectrum,model)

        # Perform F-test
        if FT.FTest(spectrum,base_model,model):
            accepted.append(i)
            accepted_models.append(model)
            accepted_names.append(param_names)
    
    # Add new component
    spectrum.p['EmissionGroups'] = AddComplexity(spectrum.p['EmissionGroups'],accepted)
    model,param_names = BD.BuildModel(spectrum)
    
    # Source starting parameters
    for source,source_params in zip(accepted_models,accepted_names):
        model = SourceParams(model,param_names,source,source_params)
    model = SourceParams(model,param_names,base_model,base_param_names)        
    
    # Fit model
    model = FitModel(spectrum,model)

    return model,param_names

# Fit Model
def FitModel(spectrum,model):
    
    # Fit model
    fit_model = fit(model,spectrum.wav,spectrum.flux,weights=spectrum.sqrtweight,maxiter=spectrum.p['MaxIter'])
    
    return fit_model

# Fit (Bootstrapped) Model
def FitBoot(spectrum,model):

    fit_model = fit(model,spectrum.wav,spectrum.Boostrap(),weights=spectrum.sqrtweight,maxiter=spectrum.p['MaxIter'])

    return fit_model.parameters

# Carry over parameters from another model
def SourceParams(model,param_names,source,source_params):
    
    # # Iterate over param names
    # for i,source_param in enumerate(source_params):
    #     if source_param in param_names:
    #         index = param_names.index(source_param)
    #         model.parameters[index] = source.parameters[i]

    # Count up number of components for a line
    numcomp = {}
    for param_name in param_names:
        if 'Flux' in param_name:
            line = param_name.split('-')[-2]
            if line not in numcomp.keys():
                numcomp[line] = 1
            else: 
                numcomp[line] += 1
    
    # Reduce flux of a line by number of components
    for i,param_name in enumerate(param_names):
        if 'Flux' in param_name:
            model.parameters[i] /= numcomp[param_name.split('-')[-2]]
            
    return model

# Add additional component to a model
def AddComplexity(EmissionGroups_old,index):

    # If multiple indices, add all of them
    if hasattr(index,'__iter__'):
        for i in index:
            EmissionGroups_old = AddComplexity(EmissionGroups_old,i)
        return EmissionGroups_old
    
    # Deepcopy
    EmissionGroups = copy.deepcopy(EmissionGroups_old)
    # Keeping track
    i = 0
    
    # Iterate in emission line dictionary 
    for group in EmissionGroups:
        for species in group['Species']:
            
            # If we have a flag
            if (species['Flag'] > 0):

                # Flag
                flag         = bin(species['Flag'])[2:]
                flag_len     = np.sum([int(bit) for bit in flag])
                
                # Check that our index is in the range
                if (index >= i) and (index < i + flag_len):
                
                    # Position in flagged bits
                    j = 0
                        
                    # Iterate over bits in flag
                    for bit in flag:
                        
                        # If we arrive at the correct index
                        if bit == '1':
                        
                            # If we've arrived at the index
                            if index == i:
                            
                                # Construct the added entry
                                entry = {
                                   'Name':species['Name'] + '-' + AC.ComponentName(j),
                                   'Lines':species['Lines'],
                                   'Flag': int('-0b1'+j*'0',2),
                                   'FlagGroups':[]
                                }

                                # Destination groupname
                                for group in EmissionGroups:
                                    if group['Name'] == species['FlagGroups'][j]:
                                        group['Species'].append(entry)

                                return EmissionGroups
                                                        
                            # Increment along the bit                            
                            i += 1
                            j += 1
                    
                i += flag_len

    