""" Build Model for Spectrum """

# Packages
import numpy as np
from scipy.optimize import LinearConstraint

# gelato supporting files
import gelato.CustomModels as CM
import gelato.AdditionalComponents as AC

# Build the model from continuum and emission
def BuildModel(spec,cont,cont_x,emis,emis_x,constraints=[]):

    return CM.CompoundModel(cont.models+emis.models,constraints=constraints),np.concatenate([cont_x,emis_x])

# Build the emission lines from the EmissionGroups and spectrum with initial guess
def BuildEmission(spectrum, EmissionGroups=None):
    
    ## Build Base Model
    model_components = []
    x0 = []

    # Check if we were passed an emission lines
    if EmissionGroups == None:
        EmissionGroups = spectrum.p['EmissionGroups']

    # Over all emission lines
    for group in EmissionGroups:
        for species in group['Species']:
            for line in species['Lines']:
                name =  group['Name']+'_'+species['Name']+'_'+str(line['Wavelength'])
                
                # If additional component
                if species['Flag'] < 0:
                    model,x = AC.AddComponent(species['Flag'],line['Wavelength'],spectrum,name)
                    x0.append(x)

                # IF original model component
                else:
                    model = CM.SpectralFeature(line['Wavelength'],spectrum,name)
                    x0.append(model.starting())

                # Add model
                model_components.append(model)
                
    ## Build Base Model

    return CM.CompoundModel(model_components),np.concatenate(x0)

# Tied Emission model parameters
def TieParams(spectrum,param_names,EmissionGroups=None):

    constraints = []

    if type(EmissionGroups) == type(None):
        EmissionGroups = spectrum.p['EmissionGroups']

    ## Tie parameters ##
    for group in EmissionGroups:
        first_group_member = True

        for species in group['Species']:
            first_species_line = True
            first_species_flux = True

            for line in species['Lines']:
                
                # Find parameter names
                name =  group['Name'] + '_' + species['Name'] + '_' + str(line['Wavelength']) + '_'
                name_Redshift = name + 'Redshift'
                name_Dispersion = name + 'Dispersion'
                name_Flux = name + 'Flux'

                ## Tie Group Components
                # Check for first line in group and make tie functions
                if first_group_member:
                    first_group_member = False
                    first_group_redshift_index = param_names.index(name_Redshift)
                    first_group_dispersion_index = param_names.index(name_Dispersion)
                else:
                    # Otherwise tie redshift (check if we should)
                    if group['TieRedshift']:
                        constraints.append(np.zeros(len(param_names)))
                        constraints[-1][first_group_redshift_index] = 1
                        constraints[-1][param_names.index(name_Redshift)] = 1
                    # Otherwise tie dispersion (check if we should)
                    if group['TieDispersion']:
                        constraints.append(np.zeros(len(param_names)))
                        constraints[-1][first_group_dispersion_index] = 1
                        constraints[-1][param_names.index(name_Dispersion)] = 1

                ## Tie Species Components
                # Tie Dispersion and Redshift
                # Check for first line in species and make tie functions
                if first_species_line:
                    # Dispersion
                    first_species_line = False
                    first_species_redshift_index = param_names.index(name_Redshift)
                    first_species_dispersion_index = param_names.index(name_Dispersion)

                # Otherwise tie params
                else:
                    # Redshift
                    constraints.append(np.zeros(len(param_names)))
                    constraints[-1][first_species_redshift_index] = 1
                    constraints[-1][param_names.index(name_Redshift)] = 1
                    # Dispersion
                    constraints.append(np.zeros(len(param_names)))
                    constraints[-1][first_species_dispersion_index] = 1
                    constraints[-1][param_names.index(name_Dispersion)] = 1

                # Tie Flux
                if (first_species_flux and (line['RelStrength'] != None)):
                    first_species_flux = False
                    # Flux
                    reference_flux = line['RelStrength']
                    index_flux = param_names.index(name_Flux)
                elif (line['RelStrength'] != None):
                    constraints.append(np.zeros(len(param_names)))
                    constraints[-1][index_flux] = 1
                    constraints[-1][param_names.index(name_Flux)] = line['RelStrength']/reference_flux

    ## Tie parameters ##

    ## Create Constraints ##
    sorton = np.array([tuple(np.argwhere(c).T[0]) for c in constraints],dtype=[('a','<i8'),('b','<i8')])
    constraints = np.array(constraints)[np.argsort(sorton)]

    # Get unique constraints
    cont_list = []
    tar_idx = []
    for c in constraints:
        ref_i = np.argwhere(c)[0][0] # Reference index
        tar_i = np.argwhere(c)[1][0] # Target index
        if tar_i not in tar_idx:
            tar_idx.append(tar_i)
            cont_list.append((ref_i,tar_i,c[tar_i]))
    ## Create Constraints ##

    return cont_list