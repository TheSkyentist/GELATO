""" Fit Model for Spectrum """

# Ignore warnings
import warnings
warnings.simplefilter('ignore')

# Packages
import copy
import numpy as np
from itertools import combinations
from astropy.modeling import fitting
fit = fitting.LevMarLSQFitter()

# gelato supporting files
import gelato.BuildModel as BM
import gelato.CustomModels as CM
import gelato.ModelComparison as MC
import gelato.AdditionalComponents as AC

# Perform initial fit of continuum with F-test
def FitContinuum(spectrum):

    # SSP Continuum
    sspcontinuum = CM.SSPContinuum(spectrum)
    ssp_names = sspcontinuum.get_names()

    # Power Law Continuum
    plcontinuum = CM.PowerLawContinuum(spectrum)
    pl_names = plcontinuum.get_names()

    # Fit initial continuuum
    ssp_fit = FitModel(spectrum,sspcontinuum,np.invert(spectrum.emission_region))

    # Try power law fit
    continuum = BM.TieContinuum(ssp_fit+plcontinuum,ssp_names+pl_names)
    continuum[1].parameters[0] *= 1/(continuum[0].parameters.size-1)
    cont_fit = FitModel(spectrum,continuum,np.invert(spectrum.emission_region))

    # Perform F-test
    if MC.FTest(spectrum,ssp_fit,cont_fit,np.invert(spectrum.emission_region)):
        
        # Fix redshift
        sspcontinuum.fix_params()
        sspcontinuum.set_region(np.ones(spectrum.wav.shape,dtype=bool))

        # Return w/ PL component
        cont = sspcontinuum+plcontinuum
        cont.parameters = cont_fit.parameters
        return cont,ssp_names+pl_names
    
    # Fix redshift
    ssp_fit.fix_params()
    ssp_fit.set_region(np.ones(spectrum.wav.shape,dtype=bool))

    return ssp_fit,ssp_names

# Construct Full Model with F-tests for additional parameters
def FitComponents(spectrum,continuum,cont_pnames,emission,emiss_pnames):

    # Base Model
    base_model = BM.BuildModel(continuum,emission)
    base_model = BM.TieParams(spectrum,base_model,cont_pnames+emiss_pnames)
    base_model = FitModel(spectrum,base_model)

    # Find number of flags
    flags = 0
    for group in spectrum.p['EmissionGroups']:
        for species in group['Species']:
            flagbits = bin(species['Flag'])[2:]
            flags += len(flagbits.replace('0',''))

    ## Use F-test for additional component selection
    # Keep track of accepted flags
    accepted = []
    # Iterate over additional components
    for i in range(flags):
        
        # Add new component
        EmissionGroups = AddComplexity(spectrum.p['EmissionGroups'],i)
        emission,emiss_pnames = BM.BuildEmission(spectrum,EmissionGroups)

        # Split FLux
        model = BM.BuildModel(continuum,emission)
        model = SplitFlux(model,cont_pnames+emiss_pnames)
        
        # Fit model
        model = BM.TieParams(spectrum,model,cont_pnames+emiss_pnames,EmissionGroups)
        model = FitModel(spectrum,model)

        # Perform F-test
        if MC.FTest(spectrum,base_model,model):
            accepted.append(i)

    ## Check all combinations of accepted components with AICs
    # All combinations
    combs = sum([list(combinations(accepted,i+1)) for i in range(len(accepted))],[])

    # Initialize AIC list
    AICs = np.zeros(len(combs))

    # Iterate over all combinations and record AICs
    for i,c in enumerate(combs):

        # Add new components
        EmissionGroups = AddComplexity(spectrum.p['EmissionGroups'],c)
        emission,emiss_pnames = BM.BuildEmission(spectrum,EmissionGroups)

        # Split Flux
        model = BM.BuildModel(continuum,emission)
        model = SplitFlux(model,cont_pnames+emiss_pnames)

        # Fit model
        model = BM.TieParams(spectrum,model,cont_pnames+emiss_pnames,EmissionGroups)
        model = FitModel(spectrum,model)

        # Calcualte AIC
        AICs[i] = MC.AIC(model,spectrum)

    # Use min AIC
    if combs != []:
        accepted = combs[np.argmin(AICs)]
    else: 
        accepted = []

    # Construct Final Model
    EmissionGroups = AddComplexity(spectrum.p['EmissionGroups'],accepted)
    emission,emiss_pnames = BM.BuildEmission(spectrum,EmissionGroups)

    # Split Flux
    model = BM.BuildModel(continuum,emission)
    model = SplitFlux(model,cont_pnames+emiss_pnames)
    # Fit model
    model = BM.TieParams(spectrum,model,cont_pnames+emiss_pnames,EmissionGroups)
    model = FitModel(spectrum,model)

    return model,cont_pnames+emiss_pnames

# Fit Model, must specify region
def FitModel(spectrum,model,region = None):

    if type(region) == type(None):
        region = np.ones(spectrum.wav.shape,dtype=bool)

    # Fit model
    fit_model = fit(model,spectrum.wav[region],spectrum.flux[region],weights=spectrum.sqrtweight[region],maxiter=spectrum.p['MaxIter'])
    
    return fit_model

# Fit (Bootstrapped) Model
def FitBoot(spectrum,model,i,N):

    # Loading bar params
    if ((spectrum.p['NProcess'] == 1) and spectrum.p['Verbose']):
        p = int(100*i/spectrum.p['NBoot']) # Percentage
        l = int(N*i/(spectrum.p['NBoot'])) # Length of bar
        if p == 0:
            print('Progress: |'+N*'-'+'|   0%',end='\r')
        elif p < 10:
            print('Progress: |'+l*'#'+(N-l)*'-'+'|   '+str(p)+'%',end='\r')
        else:
            print('Progress: |'+l*'#'+(N-l)*'-'+'|  '+str(p)+'%',end='\r')

    # Fit model
    fit_model = fit(model,spectrum.wav,spectrum.Boostrap(),weights=spectrum.sqrtweight,maxiter=spectrum.p['MaxIter'])

    return np.concatenate([fit_model.parameters,[MC.rChi2(spectrum,fit_model)]])

# Split flux between emission lines
def SplitFlux(model,param_names):

    # Count up number of components for a line
    numcomp = {}
    for param_name in param_names:
        if 'Flux' in param_name:
            if model.parameters[param_names.index(param_name)] >= 0:
                line = param_name.split('-')[-2]
                if line not in numcomp.keys():
                    numcomp[line] = 1
                else: 
                    numcomp[line] += 1
    
    # Reduce flux of a line by number of components
    for i,param_name in enumerate(param_names):
        if 'Flux' in param_name:
            n = numcomp[param_name.split('-')[-2]]
            if n > 1:
                parameters = model.parameters
                parameters[i] = model.parameters[i]/n
                model.parameters = parameters
            
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
                flag         = bin(species['Flag'])[2:][::-1]
                flag_len     = np.sum([int(bit) for bit in flag])
                
                # Check that our index is in the range
                if (index >= i) and (index < i + flag_len):
                
                    # Position in flagged bits
                    j = 0
                        
                    # Iterate over bits in flag
                    for k,bit in enumerate(flag):
                        
                        # If we arrive at the correct index
                        if bit == '1':
                        
                            # If we've arrived at the index
                            if index == i:
                            
                                # Construct the added entry
                                entry = {
                                   'Name':species['Name'] + '-' + AC.ComponentName(k),
                                   'Lines':species['Lines'],
                                   'Flag': int('-0b1'+k*'0',2),
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