""" Fit Model for Spectrum """

# # Ignore warnings
# import warnings
# warnings.simplefilter('ignore')

# Packages
import copy
import numpy as np
from itertools import combinations
from scipy.optimize import least_squares

# gelato supporting files
import gelato.BuildModel as BM
import gelato.CustomModels as CM
import gelato.ModelComparison as MC
import gelato.AdditionalComponents as AC

# Perform initial fit of continuum with F-test
def FitContinuum(spectrum):

    # Continuum region
    region = np.invert(spectrum.emission_region)
    args = (spectrum.wav[region],spectrum.flux[region],spectrum.isig[region])

    # SSP Continuum        
    ssp = CM.CompoundModel([CM.SSPContinuumFree(spectrum)])

    # Fit initial continuuum with free redshift
    x0 = ssp.starting()
    sspfit = FitModel(ssp,x0,args).x

    # SSP+PL Continuum
    pl = CM.PowerLawContinuum(spectrum,nssps=ssp.nparams()-1)
    ssppl = CM.CompoundModel([ssp.models[0],pl])

    # Starting values
    x0 = ssppl.starting()
    x0[:len(sspfit)] = sspfit

    # Fit initial continuuum with free redshift
    sspplfit = FitModel(ssppl,x0,args).x

    # Perform F-test
    if MC.FTest(ssp,sspfit,ssppl,sspplfit,spectrum,args):

        # Final Redshift (Round for Stability)
        z = sspplfit[0]

        # Get fixed redshift compound model
        sspfixed = CM.SSPContinuumFixed(z,spectrum)
        cont = CM.CompoundModel([sspfixed,pl])

        # Get starting values
        x0 = sspplfit[1:]

        # Return w/ PL component
        return cont,x0

    # Final Redshift (Round for Stability)
    z = sspfit[0]

    # Fixed Redshift compound model
    sspfixed = CM.CompoundModel([CM.SSPContinuumFixed(z,spectrum)])
    sspfixed.starting()
    
    return sspfixed,sspfit[1:]

# Construct Full Model with F-tests for additional parameters
def FitComponents(spectrum,cont,cont_x,emis,emis_x):

    # Fit region
    args = (spectrum.wav,spectrum.flux,spectrum.isig)

    # Base Model
    constraints = BM.TieParams(spectrum,cont.get_names()+emis.get_names())
    base_model,x0 = BM.BuildModel(spectrum,cont,cont_x,emis,emis_x,constraints)

    # Initial fit
    x0 = base_model.constrain(x0) # Limit to true parameters
    base_fit = FitModel(base_model,x0,args,jac=base_model.jacobian).x

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
        emis,emis_x = BM.BuildEmission(spectrum,EmissionGroups)

        # Create New Model
        constraints = BM.TieParams(spectrum,cont.get_names()+emis.get_names(),EmissionGroups)
        model,x0 = BM.BuildModel(spectrum,cont,cont_x,emis,emis_x,constraints)

        # Inital guess, split flux amongst new lines
        x0 = SplitFlux(model,x0)
        x0 = model.constrain(x0) # Limit to true parameters

        # Fit Model
        fit = FitModel(model,x0,args,jac=model.jacobian)

        # Perform F-test
        if not MC.FTest(base_model,base_fit,model,fit.x,spectrum,args):
            continue

        # Accept
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
        emis,emis_x = BM.BuildEmission(spectrum,EmissionGroups)

        # Create New Model
        constraints = BM.TieParams(spectrum,cont.get_names()+emis.get_names(),EmissionGroups)
        model,x0 = BM.BuildModel(spectrum,cont,cont_x,emis,emis_x,constraints)

        # Inital guess, split flux amongst new lines
        x0 = SplitFlux(model,x0)
        x0 = model.constrain(x0) # Limit to true parameters

        # Fit Model
        model_fit = FitModel(model,x0,args,jac=model.jacobian).x

        # Calcualte AIC
        AICs[i] = MC.AIC(model,model_fit,args)

    # Use min AIC
    if combs != []:
        accepted = combs[np.argmin(AICs)]
    else: 
        accepted = []

    # Construct Final Model
    EmissionGroups = AddComplexity(spectrum.p['EmissionGroups'],accepted)
    emis,emis_x = BM.BuildEmission(spectrum,EmissionGroups)

    constraints = BM.TieParams(spectrum,cont.get_names()+emis.get_names(),EmissionGroups)
    model,x0 = BM.BuildModel(spectrum,cont,cont_x,emis,emis_x,constraints)

    # Inital guess, split flux amongst new lines
    x0 = SplitFlux(model,x0)
    x0 = model.constrain(x0) # Limit to true parameters

    # Fit Model
    model_fit = FitModel(model,x0,args,jac=model.jacobian).x

    return model,model_fit

# Fit Model
def FitModel(model,x0,args,jac='3-point'):
    
    # fit = least_squares(fun = model.residual, jac = jac, x0 = x0, args = args, bounds = model.get_bounds(), method = 'trf', x_scale='jac')

    fit = least_squares(fun = model.residual, jac = jac, x0 = x0, args = args, bounds = model.get_bounds(), method = 'trf', x_scale='jac',max_nfev=100,tr_solver='lsmr',tr_options={'regularize':True})

    return fit

# Fit (Bootstrapped) Model
def FitBoot(model,x0,spectrum,i,N):

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
    args = spectrum.wav,spectrum.Bootstrap(i),spectrum.isig
    fit_model = FitModel(model,x0,args,jac=model.jacobian)

    return np.concatenate([fit_model.x,[np.square(model.residual(fit_model.x,spectrum.wav,spectrum.flux,spectrum.isig)).sum()]])

# Split flux between emission lines
def SplitFlux(model,x0):

    # Count up number of components for a line
    numcomp = {}
    for i,param_name in enumerate(model.get_names()):
        if 'Flux' in param_name:
            if x0[i] >= 0:
                line = param_name.split('_')[-2]
                if line not in numcomp.keys():
                    numcomp[line] = 1
                else: 
                    numcomp[line] += 1
    
    # Reduce flux of a line by number of components
    for i,param_name in enumerate(model.get_names()):
        if 'Flux' in param_name:
            n = numcomp[param_name.split('_')[-2]]
            if n > 1:
                x0[i] /= n
            
    return x0

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
                                   'Name':species['Name'] + '_' + AC.ComponentName(k),
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