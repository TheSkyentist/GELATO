""" Construct and Verify input Parameters """

import json

# Construct emission line dictionary from input JSON file.
def construct(path):

    # Open file
    file = open(path,'r')
    p = json.load(file)
    file.close()

    # Iterate through and load
    emissionLines = {}
    for group in p['EmissionLines']:
        emissionLines[group['Group']] = {}
        for species in group['Species']:
            lines = []
            for line in species['Lines']:
                lines.append((line['Wavelength'],line['RelStrength']))
            emissionLines[group['Group']][species['Name']] = (lines,species['Flag'],species['FlagGroups'])

    # Reset emission lines
    p['EmissionLines'] = emissionLines
    
    return p

# Verify if emission line dictionary is in correct format. 
def verify(emissionLines):

    # Check if Emission Lines are in dictionary
    if not type(emissionLines) == dict:
        print('Emission lines must be in dictionary.')
        return False
    
    # Check if redshift group is in a dictionary
    for z in emissionLines.keys():
        if not type(emissionLines[z]) == dict:
            print('Each redshift must be a dictionary.')
            return False
        if not type(z) == str:
            print('Redshift key must be a string.')
            return False
        
        # Check if species is in an appropriate tuple or list
        for species in emissionLines[z]:
            if not ((type(emissionLines[z][species]) == tuple) or \
                    (type(emissionLines[z][species]) == list)):
                print('Each species must be a tuple or list.')
                return False
            if not (type(species) == str):
                print('Species key must be a string.')
                return False
            if not (len(emissionLines[z][species]) == 3):
                print('Species must be length 3.')
                return False
            # Check if lines are in a list or tuple
            if not (type(emissionLines[z][species][0]) == list):
                print('Lines must be in a list.')
                return False
                
            # Check if flag is an good int
            if not (type(emissionLines[z][species][1]) == int):
                print('Lines flag must be an int.')
                return False
            if not (emissionLines[z][species][1] >= 0):
                print('Lines flag must be a positive int.')
                return False
                
            # Check if additional components are in a list 
            if not (type(emissionLines[z][species][2]) == list):
                print('Line component redshift groups must be a list.')
                return False

            # Check if each line is an appropriate tuple or list
            for line in emissionLines[z][species][0]:
                if not (type(line) == tuple) or (type(line) == list):
                    print('Each line must be a tuple or list.')
                    return False
                if not (len(line) == 2):
                    print('Each line tuple must have two elements.')
                    return False
                if not ((type(line[0]) == float) or (type(line[0]) == int)):
                    print('Line wavelength must be a float or int.')
                    return False
                if not ((type(line[1]) == float) or (type(line[1]) == int) or (line[1] == None)):
                    print('Line strength must be a float or int or None.')
                    return False
                     
        
            # Check if flag bit sum is equal to length of 
            if not (sum([int(bit) for bit in bin(emissionLines[z][species][1])[2:]]) == \
                    len(emissionLines[z][species][2])):
                print('Flag number does not match number of additional components.')
                return False
                        
    return True