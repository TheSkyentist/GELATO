""" Main Function """

# Packages
import os
import numpy as np
from datetime import datetime
from astropy.table import Table

# gelato supporting files
import gelato.Plotting as PL
import gelato.BuildModel as BM
import gelato.FittingModel as FM
import gelato.SpectrumClass as SC
import gelato.ModelComparison as MC
import gelato.EquivalentWidth as EW
import gelato.ConstructParams as CP

# Get fit parameters for spectrum
def gelato(params,path,z):

    # Load Params
    params = CP.construct(params)

    # Set seed
    # np.random.seed(params["RandomSeed"])

    # Get name of file
    name = path.split("/")[-1]

    # If it exits, skip
    if os.path.exists(params["OutFolder"]+name.replace(".fits","-results.fits")) and not params["Overwrite"]:
        if params["Verbose"]:
            print('gelato already exists:',name)
        return
    if params["Verbose"]:
        print("Making gelato for",name)

    ## Load in Spectrum ##
    if params["Verbose"]:
        print("Gathering ingredients:",name)
    spectrum = SC.Spectrum(path,z,params)
    if params["Verbose"]:
        print("Ingredients gathered:",name)

    ## Create Base Model ##
    if params["Verbose"]:
        print("Making the base:",name)
    continuum,cont_pnames = FM.FitContinuum(spectrum)
    emission,emiss_pnames = BM.BuildEmission(spectrum)
    if params["Verbose"]:
        print("Base created:",name)

    # spectrum.LimitSpectrum()

    # Check if any of the lines can be fit
    if len(spectrum.regions) > 0:

        ## Fit Additional Components ##
        if params["Verbose"]:
            print("Adding flavor:",name)
        
        model,param_names = FM.FitComponents(spectrum,continuum,cont_pnames,emission,emiss_pnames)
        param_names = param_names + ["rChi2"]
        if params["Verbose"]:
            print("Flavor added:",name)

        # Bootstrap
        if params["Verbose"]:
            print("Scooping portions (this may take a while):",name)
        parameters = np.array([FM.FitBoot(spectrum,model,i) for i in range(params["NBoot"])])
        if params["Verbose"]:
            print("Portions scooped:",name)

        ## Plotting ##
        if params["Plotting"]:
            if params["Verbose"]:
                print("Presenting gelato:",name)
            # Set model parameters to median values
            model.parameters = np.median(parameters,0)[:-1]
            PL.Plot(spectrum,model,path,param_names)
            if params["Verbose"]:
                print("gelato presented:",name)

         # Get equivalent wdiths
        if params["CalcEW"]:
            if params["Verbose"]:
                print("Measuring Texture:",name)
            parameters,param_names = EW.EquivalentWidth(spectrum,model,parameters,param_names)
            if params["Verbose"]:
                print("Measured Texture:",name)
                
    # Otherwise:
    else:
        # Just continuum
        parameters = np.hstack([continuum.parameters,[MC.rChi2(spectrum,continuum)]])
        param_names = cont_pnames + ["rChi2"]
        if params["Verbose"]:
            print("Flavour not found (no lines with spectral coverage):",name)

        ## Plotting ##
        if params["Plotting"]:
            if params["Verbose"]:
                print("Presenting gelato:",name)
            PL.PlotFig(spectrum,continuum,path,param_names)
            if params["Verbose"]:
                print("Gelato presented:",name)

    # Turn into FITS table
    parameters = Table(data=parameters,names=param_names)

    # Save results
    if params["Verbose"]:
        print("Freezing results:",name)
    parameters.write(params["OutFolder"]+name.replace(".fits","-results.fits"),overwrite=True)
    if params["Verbose"]:
        print("Results freezed:",name)

    if params["Verbose"]:
        print("GELATO finished for",name)

def header():
    print("Welcome to GELATO")
    print("Galaxy/AGN Emission Line Analysis TOol")
    print("Developed by R. E. Hviding")
    print("Started making gelato at",datetime.now())

def footer():
    print("Finished making gelato at",datetime.now())