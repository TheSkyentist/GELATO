""" Custom Models """

import numpy as np
from os import path
from lmfit import Model
from astropy.io import fits

# GELATO
# from gelato.Constants import C
C = 299792.458 # km/s (Sets the units for the returned values)

C2 = C*C
OOSQRT_2_PI = 1/np.sqrt(2*np.pi)

def func(x,Redshift,Flux,Dispersion):

    # Calculate a few things in advance
    oolcpv = 1 / (self.center * (C + Redshift)) # lam * ( C + v )
    ood = 1 / Dispersion
    exponand = ( C * ood ) * ( C * x * oolcpv - 1) 

    # Calculate Gaussian
    y = C2 * Flux * np.exp(-0.5 * exponand * exponand) * oolcpv * OOSQRT_2_PI * ood
    return y

m = Model(func)

        # Just for cleanliness
        spec = self.spec 

        # Find large region around line
        for r in spec.regions:
            if (self.center*(1+spec.z) < r[1]) and (self.center*(1+spec.z) > r[0]):
                break
        inreg = np.logical_and(spec.wav > r[0],spec.wav < r[1])

        # Find small region around line 
        lwav = self.center*(1+spec.z) # Redshift center
        lwidth = lwav*spec.p['LineRegion']/(2*C)
        inline = np.logical_and(spec.wav > lwav - lwidth,spec.wav < lwav + lwidth)

        # Find starting height and then flux
        Height = np.abs(np.max(spec.flux[inline]) - np.min(spec.flux[inreg]))
        if Height == 0: Height = np.max(spec.sigma[inline])
        Flux = Height * (1 + spec.z) * self.Dispersion * self.center / ( OOSQRT_2_PI * C )

        # Set bounds
        self.Redshift_bounds = (C*(spec.z-0.001),C*(spec.z+0.001))
        Fbound = 1.5*Flux*self.Dispersion_bounds[1]/self.Dispersion
        self.Flux_bounds = (-Fbound,Fbound)

        # Return starting value
        return np.array([spec.z*C,Flux,self.Dispersion])

print(m)

def dfunc(x,Redshift,Flux,Dispersion):

    # Calculate a few things in advance
    oocpv = 1 / (C + Redshift) # 1 / (c + v)
    oolcpv = oocpv / ( self.center ) # 1 / lam * (c + v) 
    ood = 1 / Dispersion
    exponand = ( C * ood ) * ( C * x * oolcpv - 1) 
    E2 = exponand * exponand

    d_Flux = C2 * np.exp(-0.5 * E2) * ood * oolcpv * OOSQRT_2_PI
    G = Flux * d_Flux
    d_Redshift = G * ( C2 * exponand * x * ood * oolcpv - 1 ) * oocpv
    d_Dispersion = G * ( E2 - 1 ) * ood

    return np.array([d_Redshift, d_Flux, d_Dispersion])
