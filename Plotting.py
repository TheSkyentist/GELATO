''' Plotting for Fit'''

import numpy as np
import matplotlib.pyplot as plt

# Plot figure
def Plot(outfolder,name,model,full_spectrum,regions):

    # Initialize Figure
    ncols   = len(regions)
    figname = name.split('/')[-1].split('.')[-2]
    fig     = plt.figure(figsize = (5*ncols,7))
    gs      = fig.add_gridspec(ncols=ncols,nrows=2,height_ratios=[4,1],hspace=0)

    # Sort regions list for good order
    regions.sort()

    # Background
    background = np.sum([model[i] for i in range(ncols)])

    for i,region in enumerate(regions):

        good    = np.logical_and(full_spectrum[0] < region[1],full_spectrum[0] > region[0])
        wav     = full_spectrum[0][good]   
        flux    = full_spectrum[1][good]

        # Axis to plot spectrum
        ax = fig.add_subplot(gs[0,i])

        # Plot data and model
        ax.step(wav,flux,'k')
        # Plot model components
        for m in model[ncols:]:
            ax.step(wav,background(wav) + m(wav))

        # ax.step(wav,model(wav),'r')

        # Set axis
        ax.set(xlim=region,ylim=(0,ax.get_ylim()[1]),ylabel='$F_\lambda$ [$10^{-17}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$]')
        plt.setp(ax.get_yticklabels()[0],visible = False)

        # Residual Axis
        rax = fig.add_subplot(gs[1,i],sharex = ax)
        rax.step(wav,flux - model(wav),'k')
        ymax = np.max(np.abs(rax.get_ylim()))
        rax.set(xlim=region,xlabel='Wavelength [\AA]',ylim=[-ymax,ymax],ylabel='r')


    # gs[1,0].add_subplot()
    fig.suptitle(figname)
    fig.tight_layout(rect = [0, 0, 1, 0.95])
    fig.savefig(outfolder + figname + '.pdf')
    plt.close(fig)