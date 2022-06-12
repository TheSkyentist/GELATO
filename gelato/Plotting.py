#! /usr/bin/env python

""" Plotting for Fit """

# Ignore warnings
import warnings
warnings.simplefilter('ignore')

# Packages
import numpy as np
from os import path
from astropy.io import fits
from matplotlib import pyplot
from scipy.optimize import minimize

# GELATO
import gelato.Utility as U
import gelato.CustomModels as CM
import gelato.SpectrumClass as SC

# Colorblind friendly colors
colors = np.array([(0,146,146),(182,109,255),(255,182,219),(109,182,255),(146,0,0),(36,255,36),(219,109,0)])/255

# Plot all Figures
def Plot(spectrum,model,parameters,fpath):

    for i in range(3): PlotFig(spectrum,model,parameters,fpath,plottype=i)

# Plot figure
def PlotFig(spectrum,model,parameters,fpath,plottype=0):

    # Calculate Medians
    medians = np.median(parameters,0)

    # Make figure name
    figname = U.fileName(path.split(fpath)[-1])+'-'
    if plottype == 0:
        figname += 'spec'
    elif plottype == 1:
        figname += 'fit'
    elif plottype == 2:
        figname += 'comp'

    # Dont overwrite
    if path.exists(path.join(spectrum.p['OutFolder'],figname+'.pdf')) and not spectrum.p['Overwrite']:
        if spectrum.p['Verbose']:
            print('GELATO already presented:',figname)
        return

    if plottype == 0:

        # Make figure
        fig = pyplot.figure(figsize=(15,7))
        gs = fig.add_gridspec(ncols=1,nrows=2,height_ratios=[4,1],hspace=0)

        # Get Spectrum
        wav     = spectrum.wav
        flux    = spectrum.flux
        isig    = spectrum.isig

        # Model prediction
        args = wav,flux,isig
        f = model.evaluate(medians,*args)

        # Add axes
        fax,rax = fig.add_subplot(gs[0,0]),fig.add_subplot(gs[1,0])

        # Plot Power Law
        if 'PowerLaw_Index' in model.get_names():
            continuum = CM.CompoundModel(model.models[1:2]).evaluate(medians[model.models[0].nparams:],*args)
            fax.step(wav,continuum,'k',ls='--',where='mid')

        # Plot model(s)
        # for p in parameters:
        #     fax.step(wav,model.evaluate(p,*args),'r',where='mid',alpha=0.5)
        fax.step(wav,f,'r',where='mid')

        # Subplot plotting
        subplotplot(plottype,fax,rax,spectrum,args,f)
        
    elif plottype > 0:

        # Make figure
        ncols   = len(spectrum.regions)
        fig = pyplot.figure(figsize = (5*ncols,7))
        gs = fig.add_gridspec(ncols=ncols,nrows=2,height_ratios=[4,1],hspace=0)

        # Continuum and Model
        args = spectrum.wav,spectrum.flux,spectrum.isig
        if 'PowerLaw_Index' in model.get_names():
            continuum = CM.CompoundModel(model.models[0:2]).evaluate(medians,*args)
        else: 
            continuum = CM.CompoundModel(model.models[0:1]).evaluate(medians,*args)
        f = model.evaluate(medians,*args)
        
        # Iterate over regions
        for i,region in enumerate(spectrum.regions):

            # Get Spectrum
            good    = np.logical_and(spectrum.wav < region[1],spectrum.wav > region[0])
            wav     = spectrum.wav[good]
            flux    = spectrum.flux[good]
            isig    = spectrum.isig[good]
            args    = wav,flux,isig

            # Add Axes
            fax,rax = fig.add_subplot(gs[0,i]),fig.add_subplot(gs[1,i])

            # Subplot plotting
            ymin = subplotplot(plottype,fax,rax,spectrum,args,f[good])

            # Plot Continuum
            if plottype == 1:
                fax.step(wav,continuum[good],ls='-',c='k',where='mid')
            # Plot components
            elif plottype == 2:
                init = 1
                if 'PowerLaw_Index' in model.get_names():
                    init = 2
                for j in range(init,len(model.models)):
                    m = model.models[j]
                    idx = model.indices[j]
                    cm = CM.CompoundModel([m]).evaluate(medians[idx:idx+m.nparams],*(wav,flux,isig))
                    fax.step(wav,ymin+cm,'--',c='gray')
                    # for p in parameters:
                    #     cm = CM.CompoundModel([m]).evaluate(p[idx:idx+m.nparams],*(wav,flux,isig))
                    #     fax.step(wav,ymin+cm,'--',c='gray',alpha=0.5)

    # Add title and save figure
    fig.suptitle(figname.replace('_','\_')+', $z='+str(np.round(spectrum.z,3))+'$',y=0.95)
    fig.tight_layout()
    fig.savefig(path.join(spectrum.p['OutFolder'],figname+'.pdf'))
    pyplot.close(fig)

# Add Linelabels
def subplotplot(pt,fax,rax,spectrum,args,f):
    
    # Unpack
    wav,flux,isig = args

    # Plot data
    fax.step(wav,flux,'gray',where='mid')

    # Plot model(s)
    # for p in parameters:
    #     fax.step(wav,model.evaluate(p,*args),'r',where='mid',alpha=0.5)
    fax.step(wav,f,'r',where='mid')
    
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

        # Log barrier constraints
        norm = 65 if pt == 0 else 15
        x0 = np.linspace(wav.min(),wav.max(),len(inds)+2)[1:-1] # Initial guess
        linelabellocs = minimize(logbarrier,x0,args=(wav,linelocs,norm),method='Nelder-Mead',options={'adaptive':True,'maxiter':len(inds)*750}).x

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
    rax.step(wav,(flux - f)*isig,'gray',where='mid')
    ymax = np.max(np.abs(rax.get_ylim()))
    rax.set(xlim=[wav.min(),wav.max()],xlabel=r'Observed Wavelength [\AA]',ylim=[-ymax,ymax])
    rax.set_ylabel('Deviation',fontsize=15)

    return ymin

# Log barrier constraints
def logbarrier(x,wav,linelocs,norm):
    
    y = np.concatenate([[wav.min()],x,[wav.max()]])
    return np.square(x-linelocs).sum()-np.log(y[1:]-y[:-1]+(wav.min()-wav.max())/norm).sum()

# Plot from results
def plotfromresults(params,fpath,z):

    if params["Verbose"]:
        print("Presenting GELATO:",path.split(fpath)[-1])

    ## Load in Spectrum ##
    spectrum = SC.Spectrum(fpath,z,params)

    # Get just the final bit of the path
    fpath = path.split(fpath)[-1]

    ## Load Results ##
    fname = path.join(params['OutFolder'],U.fileName(fpath))+'-results.fits'
    parameters = fits.getdata(fname,'PARAMS')
    pnames =  [n for n in parameters.columns.names if not (('EW' in n) or ('RAmp' in n) or ('PowerLaw_Scale' == n))][:-1]
    ps = np.array([parameters[n] for n in pnames]).T

    ## Create model ##
    # Add continuum
    ssp_names = [n[4:] for n in pnames if (('SSP_' in n) and (n != 'SSP_Redshift'))]
    models = [CM.SSPContinuumFree(spectrum,ssp_names = ssp_names)]
    if 'PowerLaw_Index' in pnames:
        models.append(CM.PowerLawContinuum(spectrum))
        models[-1].starting()

    if spectrum.regions != []:

        # Add spectral lines
        ind = sum([m.nparams for m in models]) # index where emission lines begin
        for i in range(ind,ps.shape[1],3):
            center = float(pnames[i].split('_')[-2])
            models.append(CM.SpectralFeature(center,spectrum))

        # Final model
        model = CM.CompoundModel(models)

        # Plot
        Plot(spectrum,model,ps,fpath)

    else:

        # Final Model
        model = CM.CompoundModel(models)

        # Plot
        PlotFig(spectrum,model,ps,fpath)

    if params["Verbose"]:
        print("GELATO presented:",fpath)