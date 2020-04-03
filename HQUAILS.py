""" Main Function """

# Packages
import os
import numpy as np
from datetime import datetime
from astropy.table import Table
import astropy.io.fits as pyfits

# HQUAILS supporting files
import Plotting as PL
import BuildModel as BM
import FittingModel as FM
import SpectrumClass as SC

# Get fit parameters for spectrum
def HQUAILS(params,path,z):

    name = path.split('/')[-1]
    if params['Verbose']:
        print("HQUAILS running for",name)

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(path,z,params)
    if params['Verbose']:
        print("Loaded spectrum:",name)

    ## Create Base Model ##
    model,param_names = BM.BuildModel(spectrum)
    if params['Verbose']:
        print("Base model created:",name)

    ## Fit Additional Components ##
    model,param_names = FM.FitComponents(spectrum,model,param_names)
    param_names = np.concatenate([param_names,['rChi2']])
    if params['Verbose']:
        print("Additional components added:",name)

    # Bootstrap
    if params['Verbose']:
        print("Beginning bootstrap (this may take a while):",name)
    parameters = np.array([FM.FitBoot(spectrum,model) for i in range(params['NBoot'])])
    if params['Verbose']:
        print("Bootstrap iterations finished:",name)

    ## Save Results
    if not os.path.exists(params['OutFolder']):
        os.mkdir(params['OutFolder'])

    Table(data=parameters,names=param_names).write(params['OutFolder']+name.replace('.fits','-results.fits'),overwrite=True)
    if params['Verbose']:
        print("Results saved:",name)

    ## Plotting ##
    if params['Plotting']:
        # Set model parameters to median values
        model.parameters = np.median(parameters,0)[:-1]
        PL.Plot(spectrum,model,path)
        if params['Verbose']:
            print("Figure saved:",name)

    if params['Verbose']:
        print("HQUAILS finished running on:",name)

def header():
    print("Welcome to HQUAILS")
    print("Handy QUAsar emIssion Line fitS (pronounced like 'quails')")
    print("Developed by R. E. Hviding")
    print("Started run at:",datetime.now())

def footer():
    print("Finished run at:",datetime.now())