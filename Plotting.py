#! /usr/bin/env python

""" Plotting for Fit """
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Colorblind friendly colors
colors = np.array([(0,146,146),(182,109,255),(255,182,219),(109,182,255),(146,0,0),(36,255,36),(219,109,0)])/255

# Plot figure
def Plot(spectrum,model,path):

    # Initialize Figure
    ncols   = len(spectrum.regions)
    figname = path.split('/')[-1].replace('.fits','')
    if os.path.exists(spectrum.p['OutFolder'] + figname + '.pdf') and not spectrum.p['Overwrite']:
        if spectrum.p['Verbose']:
            print('Figure Already Plotted:',figname)
        return
    if spectrum.p['Verbose']:
        print('Plotting Figure:',figname)
    fig     = plt.figure(figsize = (5*ncols,7))
    gs      = fig.add_gridspec(ncols=ncols,nrows=2,height_ratios=[4,1],hspace=0)

    # Continuum
    continuum = model[0]
    for i,region in enumerate(spectrum.regions):

        # Choose inside wavelength
        good    = np.logical_and(spectrum.wav < region[1],spectrum.wav > region[0])
        wav     = spectrum.wav[good]
        flux    = spectrum.flux[good]
        isig    = 1/spectrum.sigma[good]

        # Axis to plot spectrum
        fax = fig.add_subplot(gs[0,i])

        # Plot data
        fax.step(wav,flux,'gray')

        # Plot model
        # Are we plotting components?
        if spectrum.p['PlotComp']:
            # Plot components
            for j in range(1,model.n_submodels()):
                fax.step(wav,continuum(spectrum.wav)[good]+model[j](wav),'--',c=colors[(j-ncols) % len(colors)])
        else:
            fax.step(wav,model(spectrum.wav)[good],'r')

        # Axis set
        ylim = list(fax.get_ylim())
        ylim[0] = np.max((0,ylim[0]))
        ylim[1] += (ylim[1] - ylim[0])/8

        # Plot Line Names
        text_height = (np.max(flux) + ylim[1])/2
        for group in spectrum.p['EmissionGroups']:
            for species in group['Species']:
                if species['Flag'] >= 0:
                    for line in species['Lines']:
                        x = line['Wavelength']*(1+spectrum.z)
                        if ((x < region[1]) and (x > region[0])):
                            fax.text(x,text_height,species['Name'],rotation=90,fontsize=12,ha='center',va='center')

        fax.set(yticks=fax.get_yticks()[1:],xlim=region,xticks=[])
        fax.set(ylabel=r'$F_\lambda$ [$10^{-17}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$]',ylim=ylim)

        # Residual Axis
        rax = fig.add_subplot(gs[1,i])
        rax.step(wav,(flux - model(spectrum.wav)[good])*isig,'gray')
        ymax = np.max(np.abs(rax.get_ylim()))
        rax.set(xlim=region,xlabel=r'Wavelength [\AA]',ylim=[-ymax,ymax])
        rax.set_ylabel('Deviation',fontsize=15)

    # Add title and save figure
    fig.suptitle(figname.replace('_','\_')+', $z='+str(np.round(spectrum.z,3))+'$')
    fig.tight_layout()
    fig.savefig(spectrum.p['OutFolder'] + figname + '.pdf')
    plt.close(fig)

# Plot from results
def plotfromresults(params,path,z):

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(path,z,params)

    if spectrum.regions != []:

        ## Load Results ##
        parameters = fits.getdata(params['OutFolder']+path.split('/')[-1].replace('.fits','-results.fits'))
        median = np.array([np.median(parameters[n]) for n in parameters.columns.names if 'EW' not in n])[:-1]
        
        ## Create model ##
        model = []
        # Add continuum
        for region in spectrum.regions:
            model.append(CM.Continuum(params['ContinuumDeg'],region))
        # Add spectral lines
        ind = (params['ContinuumDeg']+1)*len(spectrum.regions) # index where emission lines begin
        for i in range(ind,median.size,3):
            center = float(parameters.columns.names[i].split('-')[-2])
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
    from astropy.io import fits
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
        if p['NProcess'] > 1: # Mutlithread
            import multiprocessing as mp
            pool = mp.Pool(processes=p['NProcess'])
            inputs = [(copy.deepcopy(p),o['File'],o['z']) for o in objects]
            pool.starmap(plotfromresults, inputs)
            pool.close()
            pool.join()
        else: # Single Thread
            for o in objects: plotfromresults(copy.deepcopy(p),o['File'],o['z'])
