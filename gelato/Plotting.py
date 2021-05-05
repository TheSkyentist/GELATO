#! /usr/bin/env python

""" Plotting for Fit """
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Colorblind friendly colors
colors = np.array([(0,146,146),(182,109,255),(255,182,219),(109,182,255),(146,0,0),(36,255,36),(219,109,0)])/255

# Plot all Figures
def Plot(spectrum,model,path,param_names):

    for i in range(3): PlotFig(spectrum,model,path,param_names,plottype=i)

# Plot figure
def PlotFig(spectrum,model,path,param_names,plottype=0):

    # Make figure name
    figname = path.split('/')[-1].replace('.fits','')+'-'
    if plottype == 0:
        figname += 'spec'
    elif plottype == 1:
        figname += 'fit'
    elif plottype == 2:
        figname += 'comp'

    # Dont overwrite
    if os.path.exists(spectrum.p['OutFolder'] + figname + '.pdf') and not spectrum.p['Overwrite']:
        if spectrum.p['Verbose']:
            print('Gelato already presented:',figname)
        return

    if plottype == 0:

        # Make figure
        fig = plt.figure(figsize=(15,7))
        gs = fig.add_gridspec(ncols=1,nrows=2,height_ratios=[4,1],hspace=0)

        # Get Spectrum
        wav     = spectrum.wav
        flux    = spectrum.flux
        isig    = 1/spectrum.sigma

        # Add axis
        fax = fig.add_subplot(gs[0,0])

         # Plot data
        fax.step(wav,flux,'gray')

        # Plot model
        fax.step(wav,model(spectrum.wav),'r')

        # Base Y axis on flux
        ymin = np.max([0,flux.min()])
        dy = flux.max() - ymin
        ylim = [ymin,ymin+1.3*dy] # Increase Axis size by 20%
        text_height = ymin+1.2*dy 

        # Get Line Names
        texts = {}
        for group in spectrum.p['EmissionGroups']:
            for species in group['Species']:
                if species['Flag'] >= 0:
                    for line in species['Lines']:
                        x = line['Wavelength']*(1+spectrum.z)
                        if ((x < wav[-1]) and (x > wav[0])):
                            texts[x] = species['Name']

        # Plot names
        for j,lam in enumerate(np.sort(list(texts.keys()))):
            # Unpack
            name = texts[lam]
            # Calculate position
            x = wav.min()+ (j+1)*(wav.max() - wav.min())/(len(texts)+1)
            # Text
            fax.text(x,text_height,name,rotation=90,fontsize=12,ha='center',va='center')
            # Plot Lines
            fax.plot([lam,lam,x,x],[ymin+dy*1.01,ymin+dy*1.055,ymin+dy*1.075,ymin+dy*1.12],'k-',lw=0.25)

        # Axis labels and limits
        fax.set(yticks=fax.get_yticks()[1:],xlim=[wav.min(),wav.max()],xticks=[])
        fax.set(ylabel=r'$F_\lambda$ [$10^{-17}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$]',ylim=ylim)

        # Residual Axis
        rax = fig.add_subplot(gs[1,0])
        rax.step(wav,(flux - model(spectrum.wav))*isig,'gray')
        ymax = np.max(np.abs(rax.get_ylim()))
        rax.set(xlim=[wav.min(),wav.max()],xlabel=r'Wavelength [\AA]',ylim=[-ymax,ymax])
        rax.set_ylabel('Deviation',fontsize=15)
        
    elif (plottype == 1) or (plottype == 2):

        # Make figure
        ncols   = len(spectrum.regions)
        fig = plt.figure(figsize = (5*ncols,7))
        gs = fig.add_gridspec(ncols=ncols,nrows=2,height_ratios=[4,1],hspace=0)

        # Continuum
        if 'PL_Continuum_Coefficient' in param_names:
            continuum = model[0:2]
        else: continuum = model[0]
        for i,region in enumerate(spectrum.regions):

            # Get Spectrum
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
            if plottype == 1:
                fax.step(wav,model(spectrum.wav)[good],'r')
            elif plottype == 2:
                init = 1
                if 'PL_Continuum_Coefficient' in param_names:
                    init = 2
                # Plot components
                for j in range(init,model.n_submodels()):
                    fax.step(wav,continuum(spectrum.wav)[good]+model[j](wav),'--',c=colors[(j-ncols) % len(colors)])
            fax.step(wav,continuum(spectrum.wav)[good],'k-')

            # Base Y axis on flux
            ymin = np.max([0,flux.min()])
            dy = flux.max() - ymin
            ylim = [ymin,ymin+1.3*dy] # Increase Axis size by 20%
            text_height = ymin+1.2*dy 

            # Get Line Names
            texts = {}
            for group in spectrum.p['EmissionGroups']:
                for species in group['Species']:
                    if species['Flag'] >= 0:
                        for line in species['Lines']:
                            x = line['Wavelength']*(1+spectrum.z)
                            if ((x < region[1]) and (x > region[0])):
                                texts[x] = species['Name']

            # Plot names
            for j,lam in enumerate(np.sort(list(texts.keys()))):
                # Unpack
                name = texts[lam]
                # Calculate position
                x = wav.min()+ (j+1)*(wav.max() - wav.min())/(len(texts)+1)
                # Text
                fax.text(x,text_height,name,rotation=90,fontsize=12,ha='center',va='center')
                # Plot Lines
                fax.plot([lam,lam,x,x],[ymin+dy*1.01,ymin+dy*1.055,ymin+dy*1.075,ymin+dy*1.12],'k-',lw=0.25)

            fax.set(yticks=fax.get_yticks()[1:],xlim=region,xticks=[])
            fax.set(ylabel=r'$F_\lambda$ [$10^{-17}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$]',ylim=ylim)

            # Residual Axis
            text_height = ymin+1.2*dy # Put labels halfway
            rax = fig.add_subplot(gs[1,i])
            rax.step(wav,(flux - model(spectrum.wav)[good])*isig,'gray')
            ymax = np.max(np.abs(rax.get_ylim()))
            rax.set(xlim=region,xlabel=r'Wavelength [\AA]',ylim=[-ymax,ymax])
            rax.set_ylabel('Deviation',fontsize=15)

    # Add title and save figure
    fig.suptitle(figname.replace('_','\_')+', $z='+str(np.round(spectrum.z,3))+'$',y=0.95)
    fig.tight_layout()
    fig.savefig(spectrum.p['OutFolder'] + figname + '.pdf')
    plt.close(fig)

# Plot from results
def plotfromresults(params,path,z):

    if params["Verbose"]:
        print("Presenting gelato:",path.split('/')[-1])

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(path,z,params)

    if spectrum.regions != []:

        ## Load Results ##
        fname = params['OutFolder']+path.split('/')[-1].replace('.fits','')+'-results.fits'
        parameters = fits.getdata(fname)
        median = np.array([np.median(parameters[n]) for n in parameters.columns.names if 'EW' not in n])[:-1]
        
        ## Create model ##
        # Add continuum
        model = CM.SSPContinuum(spectrum)
        model.set_region(np.ones(spectrum.wav.shape,dtype=bool))
        if 'PL_Continuum_Coefficient' in parameters.columns.names:
            model += CM.PowerLawContinuum(spectrum)

        # Add spectral lines
        ind = len(model.parameters) # index where emission lines begin
        for i in range(ind,median.size,3):
            center = float(parameters.columns.names[i].split('-')[-2])
            model += CM.SpectralFeature(center,spectrum)
        
        # Finish model and add parameters
        model.parameters = median

        # Plot
        Plot(spectrum,model,path,parameters.columns.names)

    if params["Verbose"]:
        print("Gelato presented:",path.split('/')[-1]) 

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
        # Load Objects
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
