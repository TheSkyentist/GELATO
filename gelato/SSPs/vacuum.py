#!/usr/bin/env python

# Import packages
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize

# Vaccuum to Air wavelengths (Ciddor 1996)
def vacuum_to_air(vacuum): # Input in aangstroms

    s = 1e8/(vacuum*vacuum)
    f = 1 + 5.792105E-2/(238.0185 - s) + 1.67917E-3/(57.362 - s)
    
    return vacuum/f

# Back out Ciddor using minimization algorithm
def air_to_vacuum(air): # Input in aangstroms

    if hasattr(air,'__iter__'):

        return np.array([air_to_vacuum(a) for a in air])

    return minimize(lambda x: np.square(vacuum_to_air(x)-air),air,method='Nelder-Mead').x[0]

# If you run this file, recreates SSP Vacuum Wavelengths
if __name__ == "__main__":

    # Open SSP
    h = fits.open('Ech1.30Zm0.40T00.0631_iPp0.00_baseFe_LIS5.0.fits')[0].header
    ssp_wav_air = (np.arange(h['NAXIS1']) - h['CRPIX1'] + 1)*h['CDELT1'] + h['CRVAL1']
    ssp_wav_vac = air_to_vacuum(ssp_wav_air) # Convert wavelengths

    # Save to file
    fits.HDUList([fits.PrimaryHDU(ssp_wav_vac)]).writeto('SSP_Wavelength_Vacuum.fits',overwrite=True)