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
