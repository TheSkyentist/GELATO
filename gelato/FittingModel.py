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
import gelato.Utility as U
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
    acceptPL = MC.FTest(ssp,sspfit,ssppl,sspplfit,spectrum,args)

    # Build model with Fixed Redshift
    z = sspplfit[0]
    models = [CM.SSPContinuumFixed(z,spectrum,region=region)]
    if acceptPL: models.append(pl)
    cont = CM.CompoundModel(models)

    # Fit Fixed Model
    if acceptPL: x0 = sspplfit[1:]
    else: x0 = sspfit[1:]
    sspfixedfit = FitModel(cont,x0,args)

    # Remove SSPs
    model_names = cont.constrain(cont.get_names())
    fitmask = np.invert(sspfixedfit.active_mask.astype(bool))
    sspmask = np.array(['SSP_' in n for n in model_names])
    if not np.logical_and(fitmask,sspmask).sum(): # If all SSPs are rejected, keep them all
        models = [CM.SSPContinuumFixed(z,spectrum)]
        if acceptPL: models.append(pl)
        cont = CM.CompoundModel(models)
        return cont,x0
    ssp_names = [n.replace('SSP_','') for n in model_names[np.logical_and(fitmask,sspmask)]]
    x0 = np.concatenate([sspfixedfit.x[np.logical_and(fitmask,sspmask)],sspfixedfit.x[np.invert(sspmask)]])

    # Build Model with Fixed Redshift and Reduced SSPs
    models = [CM.SSPContinuumFixed(z,spectrum,ssp_names=ssp_names)]
    if acceptPL: models.append(pl)
    cont = CM.CompoundModel(models)

    # Return continuum
    return cont,x0

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

        # Get pnames of new components
        model_pnames = model.get_names()
        diff = np.setdiff1d(model_pnames,base_model.get_names(),assume_unique=True)
        
        # Reject component if we hit limits
        # mask = fit.active_mask
        # if model.constrained: mask = model.expand(mask)
        # if np.logical_or.reduce([mask[model_pnames.index(d)] for d in diff]):
        #     continue

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
        fit = FitModel(model,x0,args,jac=model.jacobian)

        # Get pnames of new components
        model_pnames = model.get_names()
        diff = np.setdiff1d(model_pnames,base_model.get_names(),assume_unique=True)

        # Reject component if we hit limits
        mask = fit.active_mask
        if model.constrained: mask = model.expand(mask)
        if np.logical_or.reduce([mask[model_pnames.index(d)] for d in diff]):
            AICs[i] = np.inf
            continue

        # Calcualte AIC
        AICs[i] = MC.AIC(model,fit.x,args)

    # Use min AIC
    accepted = []
    if (combs != []) and (min(combs) != np.inf):
        accepted = combs[np.argmin(AICs)]

    # Construct Model with Complexity
    EmissionGroups = AddComplexity(spectrum.p['EmissionGroups'],accepted)
    emis,emis_x = BM.BuildEmission(spectrum,EmissionGroups)

    constraints = BM.TieParams(spectrum,cont.get_names()+emis.get_names(),EmissionGroups)
    model,x0 = BM.BuildModel(spectrum,cont,cont_x,emis,emis_x,constraints)

    # Inital guess, split flux amongst new lines
    x0 = SplitFlux(model,x0)
    x0 = model.constrain(x0) # Limit to true parameters

    # Fit Model
    model_fit = FitModel(model,x0,args,jac=model.jacobian)

    # Remove SSPs
    model_names = model.constrain(model.get_names())
    fitmask = np.invert(model_fit.active_mask.astype(bool))
    sspmask = np.array(['SSP_' in n for n in model_names])
    if not np.logical_and(fitmask,sspmask).sum(): # If all SSPs are rejected, keep them all
        return model,model_fit.x
    ssp_names = [n.replace('SSP_','') for n in model_names[np.logical_and(fitmask,sspmask)]]
    x0 = np.concatenate([model_fit.x[np.logical_and(fitmask,sspmask)],model_fit.x[np.invert(sspmask)]])
    
    # Build Continuum with Reduced SSPs
    cmodels = [CM.SSPContinuumFixed(model.models[0].redshift,spectrum,ssp_names=ssp_names)]
    if 'PowerLaw_Index' in model.get_names(): cmodels.append(cont.models[1])
    cont = CM.CompoundModel(cmodels)

    # Construct Final Model
    constraints = BM.TieParams(spectrum,cont.get_names()+emis.get_names(),EmissionGroups)
    model,_ = BM.BuildModel(spectrum,cont,np.array(cont_x),emis,emis_x,constraints)

    # Fit Model
    model_fit = FitModel(model,np.array(x0),args,jac=model.jacobian).x

    return model,model_fit

# Fit Model
def FitModel(model,x0,args,jac='3-point'):
    
    fit = least_squares(fun = model.residual, jac = jac, x0 = x0, args = args, bounds = model.get_bounds(), method = 'trf', x_scale='jac',max_nfev=100,tr_solver='lsmr',tr_options={'regularize':True})

    return fit

# Fit (Bootstrapped) Model
def FitBoot(model,x0,spectrum,i):

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