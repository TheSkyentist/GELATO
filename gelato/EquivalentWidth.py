""" Generate Rest Equivalent Widths """

# Import packages
import numpy as np
from os import path
from astropy.io import fits
from astropy.table import Table,hstack

# GELATO
import gelato.CustomModels as CM
import gelato.SpectrumClass as SC

# Constants
C = 299792.458 # km/s

# Calculate Equivalent Width
def EquivalentWidth(spectrum,model,parameters,param_names=None):

    if type(param_names) == type(None):
        param_names = model.get_names()

    # Continuum and Model
    args = spectrum.wav,spectrum.flux,spectrum.isig
    if 'PowerLaw_Coefficient' in param_names:
        continuum = CM.CompoundModel(model.models[0:2])
    else: 
        continuum = CM.CompoundModel(model.models[0:1])

    # Where continuum parameters start
    ind = continuum.nparams()

    # Get all continuum heights
    continua = np.ones((parameters.shape[0],len(spectrum.wav)))
    for i,params in enumerate(parameters):
        continua[i] = continuum.evaluate(params,*args)

    # Empty EW array
    REWs = np.zeros((parameters.shape[0],int((parameters.shape[1]-1-ind)/3)))

    # Iterate over Emission lines
    for i in range(ind,parameters.shape[1]-1,3):

        # Get line parameters
        param_names += (param_names[i].replace('_Redshift','_REW'),)

        center = float(param_names[i].split('_')[-2]) 
        opz = 1 + parameters[:,i] # 1 + z
        flux = parameters[:,i+1] # Line flux

        # Get line
        linewav = center*opz
        linewidth = linewav*spectrum.p['LineRegion']/(2*C)

        # Continuum height
        heights = np.ones(parameters.shape[0])
        for j,lwv,lwd,ctm in zip(range(parameters.shape[0]),linewav,linewidth,continua):
            # Continuum height region
            region = np.logical_and(spectrum.wav > lwv - lwd,spectrum.wav < lwv + lwd)
            if region.sum() == 0:
                heights[j] = np.nan
            else:
                heights[j] = np.median(ctm[region])
    
        # Get REW
        REWs[:,int((i-ind)/3)] = np.abs(flux/(heights*opz))

    # Return combined Results
    return np.hstack((parameters,REWs)),param_names

# Plot from results
def EWfromresults(params,fpath,z):

    if params["Verbose"]:
        print("Measuring texture:",path.split(fpath)[-1])

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(fpath,z,params)

    if spectrum.regions != []:

        # Load name and parameters
        fname = path.join(params['OutFolder'],path.split(fpath)[-1].replace('.fits','')+'-results.fits')
        parameters = Table.read(fname)
        names = parameters.colnames

        # Dont add if already has EWs and no overwrite
        for n in names:
            if 'EW' in n:
                if not params['Overwrite']:
                    print('Texture already measured:',path.split(fpath)[-1])
                    return
                else: 
                    parameters = parameters[[n for n in names if 'EW' not in n]]
                    names = parameters.colnames
                break
        
        # Put parameters into the form we need
        parameters = np.array([parameters[n] for n in names]).T

        ## Create model ##
        # Add continuum
        models = [CM.SSPContinuumFree(spectrum)]
        if 'PowerLaw_Coefficient' in names:
            models.append(CM.PowerLawContinuum(spectrum))
            models[-1].starting()

        # Add spectral lines
        ind = sum([m.nparams for m in models]) # index where emission lines begin
        for i in range(ind,len(names)-1,3):
            center = float(names[i].split('_')[-2])
            models.append(CM.SpectralFeature(center,spectrum))

        # Final model
        model = CM.CompoundModel(models)
        
        # Calculate REWs
        parameters,param_names = EquivalentWidth(spectrum,model,parameters,names)

        # Turn into FITS table
        parameters = Table(data=parameters,names=param_names)
        
    if params["Verbose"]:
        print("Texture measured:",path.split('/')[-1])