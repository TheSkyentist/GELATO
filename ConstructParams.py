""" Construct and Verify input Parameters """

import json
import sys

# Construct emission line dictionary from input JSON file.
def construct(path):

    # Open file
    file = open(path,'r')
    p = json.load(file)
    file.close()

    # Verify
    if not verify(p):
        print('Parameters file is not correct, exiting.')
        sys.exit(1)

    return p

# Verify if parameters json file is correct.
def verify(params):

    # Check if Emission Lines are in dictionary
    if not type(params) == dict:
        print('Parameters not loaded as a dictionary.')
        return False

    # Check that all parameters are specified
    for p in ['OutFolder', 'RegionWidth', 'LineDataWidth', 'BackgroundDeg', 'MaxIter', 'NBoot', 'FThresh', 'NProcess', 'Plotting', 'PlotComp', 'Concatenate', 'Verbose', 'EmissionGroups']:
        if not p in params.keys():
            print('Parameters does not contain parameter:',p)
            return False

    # Check that all parameters are specified and of correct types
    for p in params.keys():
        # Now check that parameters are of the correct types
        if p in ['OutFolder']:
            if not (type(params[p]) == str):
                print('Parameter',p,'must be a string.')
                return False
        elif p in ['RegionWidth','LineDataWidth']:
            if not ((type(params[p]) == float) or (type(params[p]) == int)):
                print('Parameter',p,'must be an int or a float.')
                return False
        elif p in ['BackgroundDeg','MaxIter','NBoot','NProcess']:
            if not ((type(params[p]) == int) and (params[p] > 0)):
                print('Parameter',p,'must be a positive int.')
                return False
        elif p in ['FThresh']:
            if not (((type(params[p]) == float) or (type(params[p]) == int)) and ((params[p] >= 0) and (params[p] <= 1))):
                print('Parameter',p,'must be an int or a float between 0 and 1 (inclusive).')
                return False
        elif p in ['Plotting', 'PlotComp','Concatenate','Verbose']:
            if not (type(params[p]) == bool):
                print('Parameter',p,'must be a boolean.')
                return False
        elif p in ['EmissionGroups']:
            if not type(params[p] == list):
                print('Parameter',p,'must be a list.')
                return False
        else:
            print('Additional parameter was specified.')
            return False

    # Check each group
    for group in params['EmissionGroups']:

        # Check group is a dict
        if not (type(group) == dict):
            print('Group must be a dict')
            return False
        
        # Check if all keys are in there
        for g in ['Name', 'TieRedshift', 'TieSigma', 'Species']:
            if not g in group.keys():
                print('A group does not contain parameter:',g)
                return False
        
        # Check the type of each key
        for g in group.keys():
            if g in ['Name']:
                if not (type(group[g]) == str):
                    print('Group parameter',g,'must be a string.')
                    return False
            elif g in ['TieRedshift','TieSigma']:
                if not (type(group[g]) == bool):
                    print('Group parameter',g,'must be a booolean.')
                    return False
            elif g in ['Species']:
                if not (type(group[g]) == list):
                    print('Group parameter',g,'must be a list.')
                    return False
            else:
                print('Additional group parameter was specified.')
                return False

    # Iterate over groups
    for group in params['EmissionGroups']:

        # Keep track of where
        g = 'In group ' + group['Name'] + ':'
        
        # Check each species
        for species in group['Species']:            
            
            # Check species is a dict
            if not (type(species) == dict):
                print(g,'Species must be a dict')
                return False

            # Check if all keys are in there
            for s in ['Name', 'Lines', 'Flag', 'FlagGroups']:
                if not s in species.keys():
                    print(g,'A species does not contain parameter:',s)
                    return False

            # Check the type of each key
            for s in species.keys():
                if s in ['Name']:
                    if not (type(species[s]) == str):
                        print(g,'Species parameter',s,'must be a string.')
                        return False
                elif s in ['Lines','FlagGroups']:
                    if not (type(species[s]) == list):
                        print(g,'Species parameter',s,'must be a list.')
                        return False
                elif s in ['Flag']:
                    if not ((type(species[s]) == int) and (species[s] >= 0)):
                        print(g,'Species parameter',s,'must be nonnegative int.')
                        return False
                else:
                    print(g,'Additional species parameter was specified.')
                    return False

            # Keep track of where
            s = 'In group ' + group['Name'] + ' and species ' + species['Name'] + ':'

            # Check if flag bit sum is equal to length of list
            if not (sum([int(bit) for bit in bin(species['Flag'])[2:]]) == len(species['FlagGroups'])):
                print(s,'Flag number does not match number of additional components.')
                return False
            
            # Check if each FlagGroup exists
            groupnames = [g['Name'] for g in params['EmissionGroups']]
            for flagroup in species['FlagGroups']:
                if not (flagroup in groupnames):
                    print(s,'FlagGroup',flagroup,'does not exist.')
                    return False

            # Iterate over lines
            for line in species['Lines']:
                
                # Check line is a dict
                if not (type(line) == dict):
                    print(s,'Line must be a dict.')

                # Check if all keys are in there
                for l in ['Wavelength', 'RelStrength']:
                    if not l in line.keys():
                        print(s,'A species does not contain parameter:',l)
                        return False   
            
                # Check the type of each key
                for l in line.keys():
                    if l in ['Wavelength']:
                        if not (((type(line[l]) == float) or (type(line[l]) == int)) and (line[l] > 0)):
                            print(s,'Line parameter',l,'must positive float or int.')
                            return False
                    elif l in ['RelStrength']:
                        if not ((( (type(line[l]) == float) or (type(line[l]) == int)) and (line[l] > 0)) or (line[l] == None)):
                            print(s,'Line parameter',l,'must positive float or int or none')
                            return False
    return True