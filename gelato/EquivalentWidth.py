""" Generate Rest Equivalent Widths """

# Import packages
import numpy as np
from os import path
from astropy.io import fits
from astropy.table import Table,hstack

# GELATO
import gelato.Utility as U
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
    if 'PowerLaw_Index' in param_names:
        continuum = CM.CompoundModel(model.models[0:2])
    else: 
        continuum = CM.CompoundModel(model.models[0:1])

    # Get all continuum heights
    continua = np.ones((len(parameters),len(spectrum.wav)))
    for i,params in enumerate(np.array([list(p) for p in parameters])):
        continua[i] = continuum.evaluate(params,*args)

    # Iterate over lines
    for l in ['_'.join(p.split('_')[:-1]) for p in param_names if 'Flux' in p]:

        # Line Parameters
        center = float(l.split('_')[-1])
        opz = 1 + parameters[l+'_Redshift']/C
        flux = parameters[l+'_Flux']

        # Get line
        linewav = center*opz
        linewidth = linewav*spectrum.p['LineRegion']/(2*C)

        # Continuum height
        heights = np.ones(len(parameters))
        for j,lwv,lwd,ctm in zip(range(len(parameters)),linewav,linewidth,continua):
            # Continuum height region
            region = np.logical_and(spectrum.wav > lwv - lwd,spectrum.wav < lwv + lwd)
            if region.sum() == 0:
                heights[j] = np.nan
            else:
                heights[j] = np.median(ctm[region])

        # Get REW
        parameters.add_column(flux/(heights*opz),index=parameters.colnames.index(l+'_RAmp')+1,name=l+'_REW')

    # Return combined Results
    return parameters

# Plot from results
def EWfromresults(params,fpath,z):

    if params["Verbose"]:
        print("Measuring texture:",path.split(fpath)[-1])

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(fpath,z,params)

    # Get just the final bit of the path
    fpath = path.split(fpath)[-1]

    if spectrum.regions != []:

        # Load name and parameters
        fname = path.join(params['OutFolder'],U.fileName(fpath))+'-results.fits'
        results = fits.open(fname)
        parameters = Table(results['PARAMS'].data)
        pnames = parameters.colnames

        # Remove PowerLaw Scale
        parameters = parameters[[n for n in pnames if not ('PowerLaw_Scale' in n)]]
        pnames = parameters.colnames

        # Dont add if already has EWs and no overwrite
        for n in pnames:
            if 'EW' in n:
                if not params['Overwrite']:
                    print('Texture already measured:',fpath)
                    return
                else: 
                    parameters = parameters[[n for n in pnames if not ('EW' in n)]]
                    pnames = parameters.colnames
                break

        ## Create model ##
        # Add continuum
        ssp_names = [n[4:] for n in pnames if (('SSP_' in n) and (n != 'SSP_Redshift'))]
        models = [CM.SSPContinuumFree(spectrum,ssp_names = ssp_names)]
        if 'PowerLaw_Index' in pnames:
            models.append(CM.PowerLawContinuum(spectrum))
            models[-1].starting()

        # Add spectral lines
        ind = sum([m.nparams for m in models]) # index where emission lines begin
        for i in range(ind,len(pnames)-1,3):
            center = float(pnames[i].split('_')[-2])
            models.append(CM.SpectralFeature(center,spectrum))

        # Final model
        model = CM.CompoundModel(models)
        
        # Calculate REWs
        parameters = fits.BinTableHDU(EquivalentWidth(spectrum,model,parameters,pnames))
        parameters.name = 'PARAMS'
        results['PARAMS'] = parameters
        results.writeto(fname,overwrite=True)
        
    if params["Verbose"]:
        print("Texture measured:",fpath)