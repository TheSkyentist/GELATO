""" Main Function """

# Packages
import os
import numpy as np
from datetime import datetime
from astropy.table import Table

# HQUAILS supporting files
import Plotting as PL
import BuildModel as BM
import FittingModel as FM
import SpectrumClass as SC
import EquivalentWidth as EW

# Get fit parameters for spectrum
def HQUAILS(params,path,z):

    # Get name of file
    name = path.split("/")[-1]

    # If it exits, skip
    if os.path.exists(params["OutFolder"]+name.replace(".fits","-results.fits")) and not params["Overwrite"]:
        if params["Verbose"]:
            print('Spectrum Results Already Exist:',name)
        return
    if params["Verbose"]:
        print("HQUAILS running for",name)

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(path,z,params)
    if params["Verbose"]:
        print("Loaded spectrum:",name)

    ## Create Base Model ##
    model,param_names = BM.BuildModel(spectrum)
    if params["Verbose"]:
        print("Base model created:",name)

    # Check if any of the lines can be fit
    if param_names != []:

        ## Fit Additional Components ##
        model,param_names = FM.FitComponents(spectrum,model,param_names)
        param_names = np.concatenate([param_names,["rChi2"]])
        if params["Verbose"]:
            print("Additional components added:",name)

        # Bootstrap
        if params["Verbose"]:
            print("Beginning bootstrap (this may take a while):",name)
        parameters = np.array([FM.FitBoot(spectrum,model) for i in range(params["NBoot"])])
        if params["Verbose"]:
            print("Bootstrap iterations finished:",name)

        ## Plotting ##
        if params["Plotting"]:
            model.parameters = np.median(parameters,0)[:-1]
            # Set model parameters to median values
            PL.Plot(spectrum,model,path)
            if params["Verbose"]:
                print("Figure saved:",name)

    # Otherwise:
    else:
        parameters = []
        if params["Verbose"]:
            print("No lines with spectral coverage:",name)

    Table(data=parameters,names=param_names).write(params["OutFolder"]+name.replace(".fits","-results.fits"),overwrite=True)
    if params["Verbose"]:
        print("Results saved:",name)

    if params["CalcEW"]:
        EW.EWfromresults(params,path,z)

    if params["Verbose"]:
        print("HQUAILS finished running on:",name)

def header():
    print("Welcome to HQUAILS")
    print("Handy QUAsar emIssion Line fitS (pronounced like 'quails')")
    print("Developed by R. E. Hviding")
    print("Started run at:",datetime.now())

def footer():
    print("Finished run at:",datetime.now())