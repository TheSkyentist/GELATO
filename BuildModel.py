""" Build Model for Spectrum """

# Packages
import numpy as np

# HQUAILS supporting files
import CustomModels as CM
import AdditionalComponents as AC
from astropy.modeling.models import Polynomial1D

# Build the model from the EmissionGroups and spectrum with initial guess
def BuildModel(spectrum, EmissionGroups=None):
    
    ## Build Base Model
    # True param names
    param_names    = []
    
    # Model
    model_components = []
        
    # Add background regions
    for i,region in enumerate(spectrum.regions):
        
        # Name of region
        name = 'Background-' + str(i) +'-'
        
        # Generate model
        background = CM.ContinuumBackground(spectrum.p['BackgroundDeg'],region)

        # Add starting parameters
        inregion = np.logical_and(spectrum.wav > region[0],spectrum.wav < region[1])
        background.c0 = np.median(spectrum.flux[inregion])
        
        # Add model
        model_components.append(background)
        
        # Collect param names
        for pname in model_components[-1].param_names:
            param_names.append(name+pname)

    # Check if we were passed an emission lines
    if EmissionGroups == None:
        EmissionGroups = spectrum.p['EmissionGroups']

    # Over all emission lines
    for group in EmissionGroups:
        for species in group['Species']:
            for line in species['Lines']:
                name =  group['Name'] + '-' + species['Name'] + '-' + str(line['Wavelength']) + '-'
                
                # If additional component
                if species['Flag'] < 0:
                    model = AC.AddComponent(species['Flag'],line['Wavelength'],spectrum)
                    
                # IF original model component
                else:
                    model = CM.SpectralFeature(center = line['Wavelength'],spectrum = spectrum)

                # Add model
                model_components.append(model)
                
                # Collect param names
                for pname in model_components[-1].param_names:
                    param_names.append(name+pname)
    ## Build Base Model

    ## Tie parameters ##
    model = TieParams(spectrum,np.sum(model_components),param_names,EmissionGroups)
    ## Tie parameters ##

    return model, param_names

# Tied all model parameters
def TieParams(spectrum,model,param_names,EmissionGroups):

    ## Tie parameters ##
    for group in EmissionGroups:
        first_group_member = True

        for species in group['Species']:
            first_species_line = True
            first_species_flux = True

            for line in species['Lines']:
                
                # Find parameter name prefix
                name =  group['Name'] + '-' + species['Name'] + '-' + str(line['Wavelength']) + '-'

                ## Tie Group Components
                # Check for first line in group and make tie functions
                if first_group_member:
                    first_group_member = False
                    TieGroupRedshift = GenTieFunc(param_names.index(name+'Redshift'))
                    TieGroupDispersion = GenTieFunc(param_names.index(name+'Dispersion'))
                else:
                    # Otherwise tie redshift (check if we should)
                    if group['TieRedshift']:
                        model.tied[model.param_names[param_names.index(name+'Redshift')]] = TieGroupRedshift
                    # Otherwise tie dispersion (check if we should)
                    if group['TieSigma']:
                        model.tied[model.param_names[param_names.index(name+'Dispersion')]] = TieGroupDispersion
                    
                ## Tie Species Components
                # Tie Dispersion and Redshift
                # Check for first line in species and make tie functions
                if first_species_line:
                    # Dispersion
                    first_species_line = False
                    TieSpeciesRedshift = GenTieFunc(param_names.index(name+'Redshift'))
                    TieSpeciesDispersion = GenTieFunc(param_names.index(name+'Dispersion'))

                # Otherwise tie params
                else:
                    # Dispersion and Redshift
                    model.tied[model.param_names[param_names.index(name+'Redshift')]] = TieSpeciesRedshift                
                    model.tied[model.param_names[param_names.index(name+'Dispersion')]] = TieSpeciesDispersion                

                # Tie Flux
                if (first_species_flux and (line['RelStrength'] != None)):
                    first_species_flux = False
                    # Flux
                    reference_flux = line['RelStrength']
                    index_flux = param_names.index(name+'Flux')
                elif (line['RelStrength'] != None):
                    height_index = param_names.index(name+'Flux')
                    TieFlux = GenTieFunc(index_flux,scale=line['RelStrength']/reference_flux)
                    model.tied[model.param_names[height_index]] = TieFlux
            
    ##Tie parameters ##

    return model

# Generate Tie Paramater Function, needed to preserve static values in function
# Note: tried using lambda function in place, but it didn't work. 
def GenTieFunc(index,scale=1):

    def TieFunc(model):
        return model.parameters[index]*scale

    return TieFunc