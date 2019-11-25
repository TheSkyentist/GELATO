''' Plotting for Fit'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Plot figure
def Plot(spectrum,model,path):

    # Initialize Figure
    ncols   = len(spectrum.regions)
    figname = path.split('/')[-1].replace('.fits','')
    fig     = plt.figure(figsize = (5*ncols,7))
    gs      = fig.add_gridspec(ncols=ncols,nrows=2,height_ratios=[4,1],hspace=0)

    # Background
    background = np.sum([model[i] for i in range(ncols)])

    for i,region in enumerate(spectrum.regions):

        # Choose inside wavelength
        good    = np.logical_and(spectrum.wav < region[1],spectrum.wav > region[0])
        wav     = spectrum.wav[good]
        flux    = spectrum.flux[good]
        sigma   = spectrum.sigma[good]

        # Axis to plot spectrum
        fax = fig.add_subplot(gs[0,i])

        # Plot data and model
        fax.step(wav,flux,'k')
        fax.step(wav,model(wav),'r')

        # Axis set
        ylim = list(fax.get_ylim())
        ylim[0] = np.max((0,ylim[0]))
        fax.set(ylabel='$F_\lambda$ [$10^{-17}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$]',ylim=ylim)
        plt.setp(fax.get_yticklabels()[0],visible = False)

        # Residual Axis
        rax = fig.add_subplot(gs[1,i],sharex = fax)
        rax.step(wav,flux - model(wav),'k')
        ymax = np.max(np.abs(rax.get_ylim()))
        rax.set(xlim=region,xlabel='Wavelength [\AA]',ylim=[-ymax,ymax])
        rax.set_ylabel('data $-$ model',fontsize=15)

    # Add title and save figure
    fig.suptitle(figname)
    fig.tight_layout(rect = [0, 0, 1, 0.96])
    fig.savefig(spectrum.p['OutFolder'] + figname + '.pdf')
    plt.close(fig)

# Plot from results
def plotfromresults(params,path,z):

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(path,z,params)

    ## Load Results ##
    parameters = pyfits.getdata(params['OutFolder']+path.split('/')[-1].replace('.fits','-results.fits'),1)
    median = np.array([np.median(parameters[n]) for n in parameters.columns.names])

    ## Create model ##
    model = []
    # Add background
    for region in spectrum.regions:
        model.append(CM.ContinuumBackground(params['BackgroundDeg'],region))
    # Add spectral lines
    ind = (params['BackgroundDeg']+1)*len(spectrum.regions) # index where emission lines begin
    for i in range(int((median.size - ind)/3)):
        center = float(parameters.columns.names[3*i+ind].split('_')[-2])
        model.append(CM.SpectralFeature(center,spectrum))
    
    # Finish model and add parameters
    model = np.sum(model)
    model.parameters = median

    # Plot
    Plot(spectrum,model,path)

# Plot from results
if __name__ == "__main__":

    # Import if we need them
    import sys
    import copy
    import argparse
    import astropy.io.fits as pyfits
    import CustomModels as CM
    import SpectrumClass as SC
    import ConstructParams as CP

    ## Parse Arguements to find Parameter File ##
    parser = argparse.ArgumentParser()
    parser.add_argument('Parameters', type=str, help='Path to parameters file')
    parser.add_argument('--ObjectList', type=str, help='Path to object list with paths to spectra and their redshifts.')
    parser.add_argument('--Spectrum', type=str, help='Path to spectrum.')
    parser.add_argument('--Redshift', type=float, help='Redshift of object')
    args = parser.parse_args()
    p = CP.construct(args.Parameters)
    ## Parse Arguements to find Parameter File ##

    ## Verify Emission Line Dictionary ##
    if not CP.verify(p['EmissionLines']):
        print('Unable to verify emission line dictionary, exiting.')
        sys.exit(1)
    ## Verify Emission Line Dictionary ##

    # Check if we are doing single or multi
    single = args.Spectrum != None and args.Redshift != None
    multi = args.ObjectList != None

    if single == multi:
        print('Specify either Object List XOR Spectrum and Redshift.')
        print('Both or neither were entered.')
    elif single: # One Plot
        plotfromresults(p, args.Spectrum, args.Redshift)
    elif multi: # Many plots
        # Load Obkects
        objects = np.genfromtxt(args.ObjectList,delimiter=',',dtype='U100,f8',names=['File','z'])
        if p['NPool'] > 1: # Mutlithread
            pools = mp.Pool(processes=p['NPool'])
            inputs = [(copy.deepcopy(p),o['File'],o['z']) for o in objects]
            pools.starmap(plotfromresults, inputs)
        else: # Single Thread
            for o in objects: plotfromresults(copy.deepcopy(p),o['File'],o['z'])
