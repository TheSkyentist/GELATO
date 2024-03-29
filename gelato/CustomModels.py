""" Custom Models """

from os import path
import numpy as np 
from astropy.io import fits

# GELATO
from gelato.Constants import C

C2 = C*C
OOSQRT_2_PI = 1/np.sqrt(2*np.pi)

class CompoundModel():

    """
    Combine component models
    """

    def __init__(self,models,constraints=[]):

        # List of Models
        self.models = models

        # Constraints
        self.constraints = constraints
        self.contindices = list(np.sort([c[1] for c in constraints]))
        self.constrained = len(constraints) > 0

        # Index where relevant parameters start
        self.indices = np.cumsum([0]+[m.nparams for m in self.models[:-1]])

    def nparams(self):

        return sum(m.nparams for m in self.models) - len(self.constraints)

    def starting(self):

        return np.concatenate([m.starting() for m in self.models])

    def residual(self,p,x,y,isig):

        if self.constrained: p = self.expand(p)
        return ((self.evaluate(p,x,y,isig) - y)*isig)

    def evaluate(self,p,x,y,isig):

        return np.sum([m.evaluate(p[i:i+m.nparams],x,y,isig) for i,m in zip(self.indices,self.models)],0)

    def jacobian(self,p,x,y,isig):

        if self.constrained: p = self.expand(p)

        jac = np.vstack([m.jacobian(p[i:i+m.nparams],x,y,isig) for i,m in zip(self.indices,self.models)])

        # Combine jacobian from terms
        for c in self.constraints: 
            jac[c[0]] += jac[c[1]]*c[2]

        return ((np.delete(jac,self.contindices,axis=0)*isig).T)

    def get_bounds(self):

        return np.delete(np.array(sum((m.get_bounds() for m in self.models),())).T,self.contindices,axis=1)

    def get_names(self):

        return sum((m.get_names() for m in self.models),())

    def constrain(self,p):

        # Limit 
        p = np.delete(p,self.contindices)

        return p

    def expand(self,p):

        # Expand
        p = np.insert(p,self.contindices-np.arange(len(self.contindices)),values=0)
        for c in self.constraints:
            p[c[1]] = p[c[0]]*c[2]

        return p

    def expand_multiple(self,p):

        # Expand
        p = np.insert(p,self.contindices-np.arange(len(self.contindices)),values=0,axis=1)
        for c in self.constraints:
            p[:,c[1]] = p[:,c[0]]*c[2]
        return p

class SpectralFeature():

    """
    1D Gaussian parametrized with Redshift (Scaled), Flux, Dispersion
    """

    # Default Dispersion Values #km/s
    Dispersion = 130
    Dispersion_bounds = (30,500)

    # Number of parameters of this model
    nparams = 3

    def __init__(self,center,spec,prefix=''):
        
        # Keep track
        self.spec = spec
        self.prefix = prefix

        # Set center of line
        self.center = center

    def starting(self):

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

    def evaluate(self,p,x,y,isig):

        # Gaussian Parametrized with the three variables

        Redshift,Flux,Dispersion = p # Unpack

        # Calculate a few things in advance
        oolcpv = 1/(self.center * (C + Redshift)) # lam * ( C + v )
        ood = 1/Dispersion
        exponand = ( C * ood ) * ( C * x * oolcpv - 1) 

        # Calculate Gaussian
        y = C2 * Flux * np.exp(-0.5 * exponand * exponand) * oolcpv * OOSQRT_2_PI * ood
        return y

    def jacobian(self,p,x,y,isig):

        # Calculate Jacobian Matrix

        Redshift,Flux,Dispersion = p # Unpack

        # Calculate a few things in advance
        oocpv = 1/(C + Redshift) # 1 / (c + v)
        oolcpv = oocpv / ( self.center ) # 1 / lam * (c + v) 
        ood = 1 / Dispersion
        exponand = ( C * ood ) * ( C * x * oolcpv - 1) 
        E2 = exponand * exponand

        d_Flux = C2 * np.exp(-0.5 * E2) * ood * oolcpv * OOSQRT_2_PI
        G = Flux * d_Flux
        d_Redshift = G * ( C2 * exponand * x * ood * oolcpv - 1 ) * oocpv
        d_Dispersion = G * ( E2 - 1 ) * ood

        return np.array([d_Redshift, d_Flux, d_Dispersion])

    def hessian(self,p,x,y,isig):

        # Calculate Hessian Matrices

        Redshift,Flux,Dispersion = p # Unpack

        # Calculate a few things in advance
        cpv = (C + Redshift) # C + v
        oolcpv = 1 / (self.center * cpv) # lam * (c + v) 

        ood = 1 / Dispersion
        exponand = ( C * ood ) * ( C * x * oolcpv - 1) 
        E2 = exponand * exponand
        E2m1 = E2 - 1
        c2xodlcpv = C2*x*ood*oolcpv

        # First Derivatives
        d_Flux = C2 * np.exp(-0.5 * E2) * ood * oolcpv * OOSQRT_2_PI
        G = Flux * d_Flux
        d_Redshift = G * ( C2 * exponand * x * ood * oolcpv - 1 ) / cpv
        d_Dispersion = G * E2m1 * ood

        # Diagonal terms
        d_Redshift2 = ( G / (cpv * cpv)) * ( 2 + c2xodlcpv * c2xodlcpv * E2m1 - 4 * c2xodlcpv * exponand )
        d_Flux2 = np.zeros(x.shape)
        d_Dispersion2 = G * ood * ood * ( 2 + E2 * ( E2 - 5 ) )

        # Cross terms
        d_Redshift_Flux = d_Redshift / Flux 
        d_Redshift_Dipersion = (G * exponand * ood / cpv)*(c2xodlcpv*(E2 - 3) - exponand)
        d_Flux_Dispersion = d_Dispersion / Flux

        return np.array([[d_Redshift2,d_Redshift_Flux,d_Redshift_Dipersion],\
                        [d_Redshift_Flux,d_Flux2,d_Flux_Dispersion],\
                        [d_Redshift_Dipersion,d_Flux_Dispersion,d_Dispersion2]])

    def get_bounds(self):

        return self.Redshift_bounds,self.Flux_bounds,self.Dispersion_bounds

    def get_names(self):

        return tuple(self.prefix+'_'+n for n in ['Redshift','Flux','Dispersion'])

class PowerLawContinuum():

    """
    Power Law Continuum
    """

    # Default Dispersion Values #km/s
    Coefficient_bounds = (0,np.inf)
    Index_bounds = (-np.inf,np.inf)

    # Number of parameters of this model
    nparams = 2

    def __init__(self,spec,nssps = 1):
        
        # Keep track
        self.spec = spec
        self.nssps = nssps

    def starting(self):

        # Average (positive) flux
        avg = np.nanmedian(self.spec.flux[self.spec.flux>0])

        Coefficient = avg/self.nssps
        Index = 1.5
        self.Center = np.nanpercentile(self.spec.wav,20)

        # Return starting value
        return np.array([Coefficient,Index])

    def evaluate(self,p,x,y,isig):

        # Power Law

        Coefficient,Index = p # Unpack
        
        return Coefficient*((x/self.Center)**(-Index))

    def jacobian(self,p,x,y,isig):

        # Calculate Jacobian Matrix

        Coefficient,Index = p # Unpack

        # Scaled Center
        X = (x/self.Center)

        # Derivatives
        d_Coefficient = X**(-Index)
        d_Index = -Coefficient*d_Coefficient*np.log(X)

        return np.array([d_Coefficient, d_Index])

    def hessian(self,p,x,y,isig):

        # Calculate Hessian Matrices

        Coefficient,Index = p # Unpack

        # Scaled Center
        ooC = 1/self.Center
        X = x*ooC

        # Calculate a few things in advance
        XpI = X**(-Index) # X to the power of -I
        lnX = np.log(X) 
        lnCX = np.log(C*X) 
        XpIlnX = XpI*lnX

        # Diagonal Terms
        d_Coeff2 = np.zeros(x.size)
        d_Index2 = Coefficient*XpIlnX*lnX
        d_Cen2 = Coefficient*Index*(Index-1)*XpI*ooC*ooC

        # Cross terms
        d_Coeff_Index = -XpIlnX

        return np.array([[d_Coeff2,d_Coeff_Index],\
                        [d_Coeff_Index,d_Index2,]])

    def get_bounds(self):

        return self.Coefficient_bounds,self.Index_bounds

    def get_names(self):

        return 'PowerLaw_Coefficient','PowerLaw_Index'


# SSP Continuum
class SSPContinuumFixed():

    """
    Continuum from SSPs
    """

    def __init__(self, redshift, spec, ssp_names = None, region = None):

        # Keep track
        self.spec = spec
        self.redshift = redshift

        # Get SSP Names
        ssp_dir = path.join(path.dirname(path.abspath(__file__)),'SSPs','')
        with open(ssp_dir+'continuum_models.txt','r') as f: 
            continuum_models = np.sort([l.strip() for l in f.readlines()])
        if ssp_names is None: self.ssp_names = continuum_models
        else:
            self.ssp_names = [c for c in continuum_models if np.logical_or.reduce([s in c for s in ssp_names])]

        # Number of parameters 
        self.nparams = len(self.ssp_names)

        # Get SSPs
        ssps = []
        for ssp_name in self.ssp_names:
            with fits.open(ssp_dir+ssp_name) as f:
                h = f[0].header
                flux = f[0].data
                ssps.append(flux)
        self.ssps = np.array(ssps)

        # Get SSP Wavelength, depends on Vacuum or Air Wavelengths
        if self.spec.p['VacuumWav']: 
            self.ssp_wav = fits.getdata(ssp_dir+'SSP_Wavelength_Vacuum.fits')
        else: 
            self.ssp_wav = (np.arange(h['NAXIS1']) - h['CRPIX1'] + 1)*h['CDELT1'] + h['CRVAL1']

        # Interpolate SSP
        if region is None:
            self.ssps = np.array([np.interp(self.spec.wav,self.ssp_wav*(1+redshift/C),f) for f in ssps])
        else: 
            self.ssps = np.array([np.interp(self.spec.wav[region],self.ssp_wav*(1+redshift/C),f) for f in ssps])

        # Set bounds
        self.bounds = tuple((0,np.inf) for i in range(self.nparams))

    def starting(self):

        # Set initial parameters
        x = self.spec.wav*(1+self.spec.z)
        region = np.logical_and(x > self.spec.wav.min(),x < self.spec.wav.max())
        meds = np.array([np.median(s[np.logical_and(region,s>0)]) for s in self.ssps])

        # Average (positive) flux
        avg = np.nanmedian(self.spec.flux[self.spec.flux>0])

        return avg/(len(self.ssp_names)*meds)

    def evaluate(self,p,x,y,isig):

        # Combine SSPs

        return np.dot(p.T,self.ssps)

    def jacobian(self,p,x,y,isig):

        # Calculate Jacobian Matrix

        return self.ssps
    
    def hessian(self,p,x,y,isig):

        # Calculate Hessian Matrices

        return np.zeros((self.nparams,self.nparams,spec.wav.size))

    def get_bounds(self):

        return self.bounds

    def get_names(self):

        return tuple('SSP_'+x.replace('.fits','').replace('.gz','') for x in self.ssp_names)

# SSP Continuum
class SSPContinuumFree():

    """
    Continuum from SSPs
    """

    def __init__(self, spec, ssp_names = None):

        # Keep track
        self.spec = spec
        
        # List SSPs
        ssp_dir = path.join(path.dirname(path.abspath(__file__)),'SSPs','')
        with open(ssp_dir+'continuum_models.txt','r') as f: 
            continuum_models = np.sort([l.strip() for l in f.readlines()])
        if ssp_names is None: self.ssp_names = continuum_models
        else:
            self.ssp_names = [c for c in continuum_models if np.logical_or.reduce([s in c for s in ssp_names])]

        # Number of parameters 
        self.nparams = len(self.ssp_names)+1

        # Get SSPs
        ssps = []
        for ssp_name in self.ssp_names:
            with fits.open(ssp_dir+ssp_name) as f:
                h = f[0].header
                flux = f[0].data
                ssps.append(flux)
        self.ssps = np.array(ssps)

        # Get SSP Wavelength, depends on Vacuum or Air Wavelengths
        if self.spec.p['VacuumWav']: 
            self.ssp_wav = fits.getdata(ssp_dir+'SSP_Wavelength_Vacuum.fits')
        else: 
            self.ssp_wav = (np.arange(h['NAXIS1']) - h['CRPIX1'] + 1)*h['CDELT1'] + h['CRVAL1']

        # Set bounds
        self.bounds = ((C*(spec.z-0.001),C*(spec.z+0.001)),) + tuple((0,np.inf) for i in range(self.nparams-1))

    def starting(self):

        # Set initial parameters
        x = self.ssp_wav*(1+self.spec.z)
        region = np.logical_and(x > self.spec.wav.min(),x < self.spec.wav.max())
        meds = np.array([np.median(s[np.logical_and(region,s>0)]) for s in self.ssps])

        # Average (positive) flux
        avg = np.nanmedian(self.spec.flux[self.spec.flux>0])

        return np.append([self.spec.z*C],avg/(len(self.ssp_names)*meds))

    def evaluate(self,p,x,y,isig):

        # Combine SSPs

        # Get redshift and coefficients
        z = p[0]/C
        coeffs = np.array(p[1:])
        ssps = self.ssps

        ssps = np.array([np.interp(x,self.ssp_wav*(1+z),s) for s in ssps])

        return np.dot(coeffs.T,ssps)

    def get_bounds(self):

        return self.bounds

    def get_names(self):

        return ('SSP_Redshift',) + tuple('SSP_'+x.replace('.fits','').replace('.gz','') for x in self.ssp_names)