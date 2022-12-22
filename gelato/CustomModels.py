""" Custom Models """

import numpy as np
from os import path
from astropy.io import fits
from lmfit import Model,Parameters

# GELATO
# from gelato.Constants import C
C = 299792.458 # km/s (Sets the units for the returned values)

C2 = C*C
OOSQRT_2_PI = 1/np.sqrt(2*np.pi)

"""
1D Gaussian parametrized with Redshift (Scaled), Flux, Dispersion
"""

def SpectralFeature(center,spec,prefix=''):

    # Gaussian Function
    def func(x,Redshift,Flux,Dispersion):

        # Calculate a few things in advance
        oolcpv = 1 / (center * (C + Redshift)) # lam * ( C + v )
        ood = 1 / Dispersion
        exponand = ( C * ood ) * ( C * x * oolcpv - 1) 

        # Calculate Gaussian
        y = C2 * Flux * np.exp(-0.5 * exponand * exponand) * oolcpv * OOSQRT_2_PI * ood
        return y

    # Create Model
    m = Model(func,prefix=prefix)

    # Find large region around line
    for r in spec.regions:
        if (center*(1+spec.z) < r[1]) and (center*(1+spec.z) > r[0]):
            break
    inreg = np.logical_and(spec.wav > r[0],spec.wav < r[1])

    # Find small region around line 
    lwav = center*(1+spec.z) # Redshift center
    lwidth = lwav*spec.p['LineRegion']/(2*C)
    inline = np.logical_and(spec.wav > lwav - lwidth,spec.wav < lwav + lwidth)

    # Find starting height and then flux
    height = np.abs(np.max(spec.flux[inline]) - np.min(spec.flux[inreg]))
    if Height == 0: height = np.max(spec.sigma[inline])
    flux = height * (1 + spec.z) * self.Dispersion * center / ( OOSQRT_2_PI * C )
    # fbound = 1.5*flux*self.Dispersion_bounds[1]/self.Dispersion

    # Set bounds
    self.Redshift_bounds = (C*(spec.z-0.001),C*(spec.z+0.001))
    self.Flux_bounds = (-Fbound,Fbound)

    # Parameter Hints
    m.set_param_hint('Redshift',value=spec.z*C,min=30,max=500)
    m.set_param_hint('Flux',value=flux,min=0,max=2)
    m.set_param_hint('Dispersion',value=130,min=30,max=500)

    # Gaussian Derivative
    def dfunc(x,Redshift,Flux,Dispersion):

        # Calculate a few things in advance
        oocpv = 1 / (C + Redshift) # 1 / (c + v)
        oolcpv = oocpv / center # 1 / lam * (c + v) 
        ood = 1 / Dispersion
        exponand = ( C * ood ) * ( C * x * oolcpv - 1) 
        E2 = exponand * exponand

        d_Flux = C2 * np.exp(-0.5 * E2) * ood * oolcpv * OOSQRT_2_PI
        G = Flux * d_Flux
        d_Redshift = G * ( C2 * exponand * x * ood * oolcpv - 1 ) * oocpv
        d_Dispersion = G * ( E2 - 1 ) * ood

        return np.array([d_Redshift, d_Flux, d_Dispersion])

    return m,dfunc

"""
Power Law Continuum
"""

def PowerLawContinuum(spec,nssps):

    # Power Law Function
    def func(self,x,Coefficient,Index,Scale):
        
        return Coefficient*((x/Scale)**(-Index))

    # Average (positive) flux
    avg = np.nanmedian(spec.flux[spec.flux>0])/nssps
    scale = np.nanpercentile(spec.wav,20)

    # Create Model
    m = Model(func,prefix='PowerLaw_')

    # Parameter Hints
    m.set_param_hint('Coefficient',value=avg,min=0)
    m.set_param_hint('Index',value=1.5)
    m.set_param_hint('Scale',value=scale,vary=False)

    def dfunc(self,x,Coefficient,Index,Scale):

        # Scaled Center
        X = (x/Scale)

        # Derivatives
        d_Coefficient = X**(-Index)
        d_Index = -Coefficient*d_Coefficient*np.log(X)
        d_Scale = Coefficient*Index*d_Coefficient/Scale

        return np.array([d_Coefficient, d_Index, d_Scale])

    return m,dfunc

# SSP Continuum w/ Free Redshifts
def SSPContinuumFree(spec,ssp_names=None):

    """
    Continuum from SSPs
    """

    # List SSPs
    ssp_dir = path.join(path.dirname(path.abspath(__file__)),'SSPs','')
    with open(ssp_dir+'continuum_models.txt','r') as f: 
        continuum_models = np.sort([l.strip() for l in f.readlines()])
    if ssp_names is None: ssp_names = continuum_models
    else:
        ssp_names = [c for c in continuum_models if np.logical_or.reduce([s in c for s in ssp_names])]

    # Get SSPs
    ssps = []
    for ssp_name in ssp_names:
        with fits.open(ssp_dir+ssp_name) as f:
            h = f[0].header
            flux = f[0].data
            ssps.append(flux)
    ssps = np.array(ssps)

    # Get SSP Wavelength, depends on Vacuum or Air Wavelengths
    if True: #  spec.p['VacuumWav']
        ssp_wav = fits.getdata(ssp_dir+'SSP_Wavelength_Vacuum.fits')
    else: 
        ssp_wav = (np.arange(h['NAXIS1']) - h['CRPIX1'] + 1)*h['CDELT1'] + h['CRVAL1']

    # SSP Interpolation and summation
    def func(x,Redshift,**params):

        ssp_interp = np.array([np.interp(x,ssp_wav*(1+Redshift/C),s) for s in ssps])
        
        return np.dot(np.array([params[k] for k in params.keys()]),ssp_interp)
    
    # Create Model
    model = Model(func)

    # Get median flux of each SSP
    ssp_wav_z = ssp_wav*(1+spec.z)
    region = np.logical_and(ssp_wav_z > spec.wav.min(),ssp_wav_z < spec.wav.max())
    meds = np.array([np.median(s[np.logical_and(region,s>0)]) for s in ssps])

    # Get starting coefficients
    coeffs = np.nanmedian(spec.flux[spec.flux>0])/(len(ssp_names)*meds)

    # Make Parameters
    params = Parameters()
    params.add('Redshift',value=spec.z*C,min=C*(spec.z-0.001),max=C*(spec.z+0.001))
    for i,c in enumerate(coeffs): params.add(f'SSP{i}',value=c,min=0)

    return model,params