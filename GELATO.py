""" Main Function """

# Packages
import os
import numpy as np
from datetime import datetime
from astropy.table import Table

# GELATO supporting files
import Plotting as PL
import BuildModel as BM
import FittingModel as FM
import SpectrumClass as SC
import EquivalentWidth as EW

# Get fit parameters for spectrum
def GELATO(params,path,z):

    # Get name of file
    name = path.split("/")[-1]

    # If it exits, skip
    if os.path.exists(params["OutFolder"]+name.replace(".fits","-results.fits")) and not params["Overwrite"]:
        if params["Verbose"]:
            print('GELATO already exists:',name)
        return
    if params["Verbose"]:
        print("Making GELATO for",name)

    ## Load in Spectrum ##
    if params["Verbose"]:
        print("Gathering ingredients:",name)
    spectrum = SC.Spectrum(path,z,params)
    if params["Verbose"]:
        print("Ingredients gathered:",name)

    ## Create Base Model ##
    if params["Verbose"]:
        print("Making the base:",name)
    model,param_names = BM.BuildModel(spectrum)
    if params["Verbose"]:
        print("Base created:",name)

    # Check if any of the lines can be fit
    if param_names != []:

        ## Fit Additional Components ##
        if params["Verbose"]:
            print("Adding flavor:",name)
        model,param_names = FM.FitComponents(spectrum,model,param_names)
        param_names = np.concatenate([param_names,["rChi2"]])
        if params["Verbose"]:
            print("Flavor added:",name)

        # Bootstrap
        if params["Verbose"]:
            print("Scooping portions:",name)
        parameters = np.array([FM.FitBoot(spectrum,model) for i in range(params["NBoot"])])
        if params["Verbose"]:
            print("Portions scooped:",name)

        ## Plotting ##
        if params["Plotting"]:
            if params["Verbose"]:
                print("Presenting GELATO:",name)
            model.parameters = np.median(parameters,0)[:-1]
            # Set model parameters to median values
            PL.Plot(spectrum,model,path)
            if params["Verbose"]:
                print("GELATO presented:",name)

    # Otherwise:
    else:
        parameters = []
        if params["Verbose"]:
            print("Flavour not found (no lines with spectral coverage):",name)

    if params["Verbose"]:
        print("Freezing results:",name)
    Table(data=parameters,names=param_names).write(params["OutFolder"]+name.replace(".fits","-results.fits"),overwrite=True)
    if params["Verbose"]:
        print("Results freezed:",name)

    if params["CalcEW"]:
        EW.EWfromresults(params,path,z)

    if params["Verbose"]:
        print("GELATO finished for",name)

def header():
    print("Welcome to GELATO")
    print("Galaxy/AGN Emission Line Analysis TOol")
    print("Developed by R. E. Hviding")
    print("Started making GELATO at",datetime.now())

def footer():
    print("Finished making GELATO at",datetime.now())