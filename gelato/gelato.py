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
    if type(params) == str: params = CP.construct(params)

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

    ## Fit Continuum ##
    if params["Verbose"]:
        print("Making the base:",name)
    cont,cont_x = FM.FitContinuum(spectrum)
    if params["Verbose"]:
        print("Base created:",name)

    # Check if any of the lines can be fit
    if len(spectrum.regions) > 0:

        ## Fit Additional Components ##
        if params["Verbose"]:
            print("Adding flavor:",name)
        
        # Build emission line model
        emis,emis_x = BM.BuildEmission(spectrum)
        model,model_fit = FM.FitComponents(spectrum,cont,cont_x,emis,emis_x)
        if params["Verbose"]:
            print("Flavor added:",name)

        # Bootstrap
        if params["Verbose"]:
            print("Scooping portions (this may take a while):",name)
        N = 40 # Max length of progress bar
        parameters = np.array([FM.FitBoot(model,model_fit,spectrum,i,N=N) for i in range(params["NBoot"])]).astype(params['Precision'])
        if ((params['NProcess'] == 1) and params['Verbose']):
            print('Progress: |'+N*'#'+'| 100%')
        if params["Verbose"]:
            print("Portions scooped:",name)

        # Expand parameters and names
        parameters = model.expand_multiple(parameters)
        param_names = model.get_names()+('rChi2',)

        ## Plotting ##
        if params["Plotting"]:
            if params["Verbose"]:
                print("Presenting gelato:",name)
            # Set model parameters to median values
            medians = np.median(parameters,0)[:-1]
            PL.Plot(spectrum,model,medians,path)
            if params["Verbose"]:
                print("gelato presented:",name)

        # Remove redshift scaling from parameters
        ind = cont.nparams()
        for i,j in zip(range(ind,parameters.shape[1]-1,3),range(len(cont.models),len(model.models))):
            parameters[:,i] /= model.models[j].zscale

        ## Equivalent Widths ##
        if params["CalcEW"]:
            if params["Verbose"]:
                print("Measuring texture:",name)
            parameters,param_names = EW.EquivalentWidth(spectrum,model,parameters,param_names)
            if params["Verbose"]:
                print("Measured texture:",name)

    # Otherwise:
    else:

        if params["Verbose"]:
            print("Flavour not found (no lines with spectral coverage):",name)

        # Just continuum
        model = cont

        # Get Chi2
        region = np.invert(spectrum.emission_region)
        args = (spectrum.wav,spectrum.flux,spectrum.isig)
        residuals = model.residual(cont_x,*args)[region]
        parameters = np.array([np.hstack([cont_x,[np.square(residuals).sum()]])])
        param_names = model.get_names()+('rChi2',)

        ## Plotting ##
        if params["Plotting"]:
            if params["Verbose"]:
                print("Presenting gelato:",name)
                medians = np.median(parameters,0)[:-1]
            PL.PlotFig(spectrum,model,medians,path)
            if params["Verbose"]:
                print("Gelato presented:",name)

    # Add in continuum redshifts
    parameters = np.hstack([np.ones((len(parameters),1),dtype=params['Precision'])*model.models[0].redshift,parameters])
    param_names = ('SSP_Redshift',) + param_names

    # Turn Chi2 into rChi2
    parameters[:,param_names.index('rChi2')] /= len(spectrum.wav) - model.nparams()

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

    # Return model (If not multiprocessing)
    if params["NProcess"] == 1: return model

def header():
    print("Welcome to GELATO")
    print("Galaxy/AGN Emission Line Analysis TOol")
    print("Developed by R. E. Hviding")
    print("Started making gelato at",datetime.now())

def footer():
    print("Finished making gelato at",datetime.now())