""" Main Function """

# Packages
import numpy as np
from os import path
from datetime import datetime
from astropy.io import fits
from astropy.table import Table

# gelato supporting files
from gelato.Constants import C
import gelato.Plotting as PL
import gelato.BuildModel as BM
import gelato.CustomModels as CM
import gelato.FittingModel as FM
import gelato.SpectrumClass as SC
import gelato.ModelComparison as MC
import gelato.EquivalentWidth as EW
import gelato.ConstructParams as CP

# Get fit parameters for spectrum
def gelato(params,spath,z):

    # Load Params
    if type(params) == str: params = CP.construct(params)

    # Get name of file
    name = path.split(spath)[-1]

    # If it exits, skip
    if path.exists(params["OutFolder"]+name.replace(".fits","-results.fits")) and not params["Overwrite"]:
        if params["Verbose"]:
            print('GELATO already exists:',name)
        return
    if params["Verbose"]:
        print("Making GELATO for",name)

    ## Load in Spectrum ##
    if params["Verbose"]:
        print("Gathering ingredients:",name)
    spectrum = SC.Spectrum(spath,z,params)
    if params["Verbose"]:
        print("Ingredients gathered:",name)

    ## Fit Continuum ##
    if params["Verbose"]:
        print("Making the base:",name)
    try: cont,cont_x = FM.FitContinuum(spectrum)
    except np.linalg.LinAlgError:
        if params["Verbose"]: print("\nGELATO failed for:",name)
        return

    if params["Verbose"]:
        print("Base created:",name)

    # Check if any of the lines can be fit
    if len(spectrum.regions) > 0:

        ## Fit Additional Components ##
        if params["Verbose"]:
            print("Adding flavor:",name)
        
        # Build emission line model
        emis,emis_x = BM.BuildEmission(spectrum)
        try: model,model_fit = FM.FitComponents(spectrum,cont,cont_x,emis,emis_x)
        except np.linalg.LinAlgError:
            if params["Verbose"]: print("\nGELATO failed for:",name)
            return
        if params["Verbose"]:
            print("Flavor added:",name)

        # Bootstrap
        if params["NBoot"] > 1:
            if params["Verbose"]:
                print("Scooping portions (this may take a while):",name)
            N = 40 # Max length of progress bar
            parameters = np.ones((params["NBoot"],len(model_fit)+1))
            for i in range(params["NBoot"]):
                try: parameters[i] = FM.FitBoot(model,model_fit,spectrum,i,N=N)
                except np.linalg.LinAlgError:
                    if params["Verbose"]: print("\nGELATO failed for:",name)
                    return
            if ((params['NProcess'] == 1) and params['Verbose']):
                print('Progress: |'+N*'#'+'| 100%')
            if params["Verbose"]:
                print("Portions scooped:",name)
        else: 
            parameters = np.ones((1,len(model_fit)+1))
            parameters[0] = np.concatenate([model_fit,[np.square(model.residual(model_fit,spectrum.wav,spectrum.flux,spectrum.isig)).sum()]])

        # Expand parameters and names
        if model.constrained: parameters = model.expand_multiple(parameters)
        param_names = model.get_names()+('rChi2',)

        ## Plotting ##
        if params["Plotting"]:
            if params["Verbose"]:
                print("Presenting GELATO:",name)
            PL.Plot(spectrum,model,parameters,spath)
            if params["Verbose"]:
                print("GELATO presented:",name)

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
                print("Presenting GELATO:",name)
            PL.PlotFig(spectrum,model,parameters,spath)
            if params["Verbose"]:
                print("GELATO presented:",name)

    ### Save Model(s)
    # Median
    median = np.nanmedian(parameters[:,:-1],0)
    # Total Model
    total_med = model.evaluate(median,spectrum.wav,spectrum.flux,spectrum.isig)
    # Start HDUL
    hdul = [fits.PrimaryHDU()]
    # SSP Continuum
    continuum_med = CM.CompoundModel(model.models[0:1]).evaluate(median,spectrum.wav,spectrum.flux,spectrum.isig)
    # Start Median Table
    medtab = [np.log10(spectrum.wav),spectrum.flux,spectrum.weight,total_med,continuum_med]
    medtabnames = ['loglam','flux','ivar','MODEL','SSP']
    # PL Continuum
    if 'PowerLaw_Coefficient' in model.get_names():
        pl_med = CM.CompoundModel(model.models[0:2]).evaluate(median,spectrum.wav,spectrum.flux,spectrum.isig) - continuum_med
        medtab.append(pl_med)
        medtabnames.append('PL')
        continuum_med = continuum_med + pl_med
    if len(spectrum.regions) > 0:
        lines_med = total_med - continuum_med
        medtab.append(lines_med)
        medtabnames.append('LINE')
    # Finish up table
    hdul.append(fits.BinTableHDU(Table(medtab,names=medtabnames)))
    hdul[-1].name = 'SUMMARY'

    # Add in continuum redshifts
    parameters = np.hstack([np.ones((len(parameters),1))*model.models[0].redshift,parameters])
    param_names = ('SSP_Redshift',) + param_names

    # Turn Chi2 into rChi2
    parameters[:,param_names.index('rChi2')] /= len(spectrum.wav) - model.nparams()

    # Turn into FITS table
    parameters = Table(data=parameters,names=param_names)

    # Add rest line amplitudes
    for l in ['_'.join(p.split('_')[:-1]) for p in param_names if 'Flux' in p]:
        center = float(l.split('_')[-1])
        ramp = parameters[l+'_Flux']*C/(parameters[l+'_Dispersion']*center*np.sqrt(2*np.pi))
        parameters.add_column(ramp,index=parameters.colnames.index(l+'_Dispersion')+1,name=l+'_RAmp')

    ## Equivalent Widths ##
    if params["CalcEW"]:
        if params["Verbose"]:
            print("Measuring texture:",name)
        parameters = EW.EquivalentWidth(spectrum,model,parameters,param_names)
        if params["Verbose"]:
            print("Measured texture:",name)

    # Save results
    if params["Verbose"]:
        print("Freezing results:",name)
    if 'PowerLaw_Coefficient' in model.get_names():
        parameters.add_column(model.models[1].Center*np.ones(len(parameters)),index=parameters.colnames.index('PowerLaw_Index')+1,name='PowerLaw_Scale')
    hdul.append(fits.BinTableHDU(parameters))
    hdul[-1].name = 'PARAMS'
    fits.HDUList(hdul).writeto(path.join(params["OutFolder"],name.replace(".fits","-results.fits")),overwrite=True)
    if params["Verbose"]:
        print("Results freezed:",name)

    if params["Verbose"]:
        print("GELATO finished for",name)

    # Return model (If not multiprocessing)
    if params["NProcess"] == 1: return model

def loadObjects(tpath):
    if tpath.endswith('.csv'):
        objects = np.atleast_1d(np.genfromtxt(tpath,delimiter=',',dtype=['U100',np.float_],names=['Path','z']))
    elif tpath.endswith('.fits'):
        objects = Table.read(tpath)
        objects.convert_bytestring_to_unicode()
        objects = np.atleast_1d(objects)
    else:
        print('Object list not .csv or .fits.')
        exit()
    return objects

def header():
    print("Welcome to GELATO")
    print("Galaxy/AGN Emission Line Analysis TOol")
    print("Developed by R. E. Hviding")
    print("Started making gelato at",datetime.now())

def footer():
    print("Finished making GELATO at",datetime.now())