#! /usr/bin/env python

""" Plotting for Fit """

# Ignore warnings
import warnings
warnings.simplefilter('ignore')

# Packages
import os
import numpy as np
from astropy.io import fits
from matplotlib import pyplot
from scipy.optimize import minimize

# GELATO
import gelato.CustomModels as CM
import gelato.SpectrumClass as SC

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
        fig = pyplot.figure(figsize=(15,7))
        gs = fig.add_gridspec(ncols=1,nrows=2,height_ratios=[4,1],hspace=0)

        # Get Spectrum
        wav     = spectrum.wav
        flux    = spectrum.flux
        isig    = 1/spectrum.sigma

        # Add axis
        fax = fig.add_subplot(gs[0,0])

         # Plot data
        fax.step(wav,flux,'gray',where='mid')

        # Plot model
        fax.step(wav,model(spectrum.wav),'r',where='mid')

        # Base Y axis on flux
        ymin = np.max([0,flux.min()])
        dy = flux.max() - ymin
        ylim = [ymin,ymin+1.3*dy] # Increase Axis size by 20%
        text_height = ymin+1.2*dy 

        # Get Line Names/Positions
        linelocs = []
        linelabels = []
        for group in spectrum.p['EmissionGroups']:
            for species in group['Species']:
                if species['Flag'] >= 0:
                    for line in species['Lines']:
                        x = line['Wavelength']*(1+spectrum.z)
                        if ((x < wav[-1]) and (x > wav[0])):
                            linelocs.append(x)
                            linelabels.append(species['Name'])

        # If we have lines to plot
        if len(linelocs) > 0:

            # Reorder line positions
            inds = np.argsort(linelocs)
            linelocs = np.array(linelocs)[inds]
            linelabels = np.array(linelabels)[inds]
            linelabellocs = minimize(lambda x: np.square(x-linelocs).sum()-np.log(x[1:]-x[:-1]+(wav.min()-wav.max())/65).sum(),np.linspace(linelocs.min(),linelocs.max(),len(inds)),method='Nelder-Mead',options={'adaptive':True,'maxiter':len(inds)*750}).x

            # Plot names
            for lineloc,linelabel,linelabelloc in zip(linelocs,linelabels,linelabellocs):
                # Text
                fax.text(linelabelloc,text_height,linelabel,rotation=90,fontsize=12,ha='center',va='center')
                # Plot Lines
                fax.plot([lineloc,lineloc,linelabelloc,linelabelloc],[ymin+dy*1.01,ymin+dy*1.055,ymin+dy*1.075,ymin+dy*1.12],ls='-',c='gray',lw=0.25)

        # Axis labels and limits
        fax.set(ylabel=r'$F_\lambda$ ['+spectrum.p['FlamUnits']+']',ylim=ylim)
        fax.set(yticks=[t for t in fax.get_yticks() if (t > ymin+0.05*dy) and (t < ylim[-1])],xlim=[wav.min(),wav.max()],xticks=[])

        # Residual Axis
        rax = fig.add_subplot(gs[1,0])
        rax.step(wav,(flux - model(spectrum.wav))*isig,'gray',where='mid')
        ymax = np.max(np.abs(rax.get_ylim()))
        rax.set(xlim=[wav.min(),wav.max()],xlabel=r'Wavelength [\AA]',ylim=[-ymax,ymax])
        rax.set_ylabel('Deviation',fontsize=15)
        
    elif (plottype == 1) or (plottype == 2):

        # Make figure
        ncols   = len(spectrum.regions)
        fig = pyplot.figure(figsize = (5*ncols,7))
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
            fax.step(wav,flux,'gray',where='mid')

            # Plot model
            # Are we plotting components?
            if plottype == 1:
                fax.step(wav,model(spectrum.wav)[good],'r',where='mid')
            elif plottype == 2:
                init = 1
                if 'PL_Continuum_Coefficient' in param_names:
                    init = 2
                # Plot components
                for j in range(init,model.n_submodels()):
                    fax.step(wav,continuum(spectrum.wav)[good]+model[j](wav),'--',c=colors[(j-ncols) % len(colors)])
            fax.step(wav,continuum(spectrum.wav)[good],ls='-',c='k',where='mid')

            # Base Y axis on flux
            ymin = np.max([0,flux.min()])
            dy = flux.max() - ymin
            ylim = [ymin,ymin+1.3*dy] # Increase Axis size by 20%
            text_height = ymin+1.2*dy 
                
            # Get Line Names/Positions
            linelocs = []
            linelabels = []
            for group in spectrum.p['EmissionGroups']:
                for species in group['Species']:
                    if species['Flag'] >= 0:
                        for line in species['Lines']:
                            x = line['Wavelength']*(1+spectrum.z)
                            if ((x < wav[-1]) and (x > wav[0])):
                                linelocs.append(x)
                                linelabels.append(species['Name'])

            # If we have lines to plot
            if len(linelocs) > 0:

                # Reorder line positions
                inds = np.argsort(linelocs)
                linelocs = np.array(linelocs)[inds]
                linelabels = np.array(linelabels)[inds]
                linelabellocs = minimize(lambda x: np.square(x-linelocs).sum()-np.log(x[1:]-x[:-1]+(wav.min()-wav.max())/15).sum(),np.linspace(wav.min(),wav.max(),len(inds)+2)[1:-1],method='Nelder-Mead',options={'adaptive':True,'maxiter':len(inds)*500}).x

                # Plot names
                for lineloc,linelabel,linelabelloc in zip(linelocs,linelabels,linelabellocs):
                    # Text
                    fax.text(linelabelloc,text_height,linelabel,rotation=90,fontsize=12,ha='center',va='center')
                    # Plot Lines
                    fax.plot([lineloc,lineloc,linelabelloc,linelabelloc],[ymin+dy*1.01,ymin+dy*1.055,ymin+dy*1.075,ymin+dy*1.12],'k-',lw=0.25)

            fax.set(ylabel=r'$F_\lambda$ ['+spectrum.p['FlamUnits']+']',ylim=ylim)
            fax.set(yticks=[t for t in fax.get_yticks() if (t > ymin+0.05*dy) and (t < ylim[-1])],xlim=[wav.min(),wav.max()],xticks=[])

            # Residual Axis
            text_height = ymin+1.2*dy # Put labels halfway
            rax = fig.add_subplot(gs[1,i])
            rax.step(wav,(flux - model(spectrum.wav)[good])*isig,'gray',where='mid')
            ymax = np.max(np.abs(rax.get_ylim()))
            rax.set(xlim=region,xlabel=r'Wavelength [\AA]',ylim=[-ymax,ymax])
            rax.set_ylabel('Deviation',fontsize=15)

    # Add title and save figure
    fig.suptitle(figname.replace('_','\_')+', $z='+str(np.round(spectrum.z,3))+'$',y=0.95)
    fig.tight_layout()
    fig.savefig(spectrum.p['OutFolder'] + figname + '.pdf')
    pyplot.close(fig)

# Plot from results
def plotfromresults(params,path,z):

    if params["Verbose"]:
        print("Presenting gelato:",path.split('/')[-1])

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(path,z,params)

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

    if spectrum.regions != []:

        # Add spectral lines
        ind = len(model.parameters) # index where emission lines begin
        for i in range(ind,median.size,3):
            center = float(parameters.columns.names[i].split('-')[-2])
            model += CM.SpectralFeature(center,spectrum)
        
        # Finish model and add parameters
        model.parameters = median

        # Plot
        Plot(spectrum,model,path,parameters.columns.names)

    else:

        # Finish model and add parameters
        model.parameters = median

        # Plot
        PlotFig(spectrum,model,path,parameters.columns.names)

    if params["Verbose"]:
        print("Gelato presented:",path.split('/')[-1])