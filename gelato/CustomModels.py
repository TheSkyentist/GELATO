""" Custom Astropy Models """

import os
import numpy as np
from spectres import spectres
from scipy.signal import fftconvolve
from astropy.io import fits
from astropy.modeling import Fittable1DModel,Parameter
from astropy.modeling.polynomial import PolynomialModel

# Constants
FLOAT_EPSILON = float(np.finfo(np.float32).tiny)
SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))
SQRT_2_PI = np.sqrt(2*np.pi)
C = 299792.458 # km/s

## Custom Astropy Modeling Models
# Spectral Feature Class (Gaussian)
class SpectralFeature(Fittable1DModel):

    """
    Very similar to the 1D Gaussian but parametrized differently. 
    """

    # Parameters
    Redshift   = Parameter(default=0)
    Flux       = Parameter(default=0,bounds=(0,None)) # Must be non-negative for emission
    Dispersion = Parameter(default=150,bounds=(60,500)) # Reasonable range

    def __init__(self, center, spectrum, Dispersion = Dispersion.default, **kwargs):
        
        # Find corresponding region
        for region in spectrum.regions:
            if (center*(1+spectrum.z) < region[1]) and (center*(1+spectrum.z) > region[0]):
                break
        inregion = np.logical_and(spectrum.wav > region[0],spectrum.wav < region[1])

        # Find region associated with emission line
        linewav = center*(1+spectrum.z)
        linewidth = linewav*spectrum.p['LineRegion']/(2*C)

        inline = np.logical_and(spectrum.wav > linewav - linewidth,spectrum.wav < linewav + linewidth)

        # Find starting height and then flux
        Height = np.abs(np.max(spectrum.flux[inline]) - np.median(spectrum.flux[inregion]))
        Flux = Height * (1 + spectrum.z) * Dispersion * center * SQRT_2_PI / C

        # Set parameters
        self.center = center
        self.domain = region
        super().__init__(Redshift = spectrum.z, Flux=Flux, Dispersion=Dispersion, **kwargs)
        self.Redshift.bounds = (spectrum.z - 0.005,spectrum.z + 0.005)
        self.Flux.bounds = (0,1.5*Flux*self.Dispersion.bounds[1]/self.Dispersion) # Set positive bounds

    @property
    def sigma(self):
        """Gaussian full Sigma at half maximum."""
        return self.Dispersion * self.center / C 
        
    @property
    def fwhm(self):
        """Gaussian full Sigma at half maximum."""
        return self.sigma() * SIGMA_TO_FWHM

    def evaluate(self, x, Redshift, Flux, Dispersion):
        """
        Gaussian1D model function.
        """

        # (1 + z) * lambda
        lam_obs = (1 + Redshift) * self.center
        exponand = (C / Dispersion) * (x / lam_obs - 1) 
        y = C * Flux * np.exp(-0.5 * exponand * exponand) / (lam_obs * Dispersion * SQRT_2_PI)    

        # # Zero outside region
        # indomain = np.logical_and(x > self.domain[0],x < self.domain[1])        
        
        # # Assign values
        # y = np.zeros(x.shape)
        # exponand = (C / Dispersion) * (x[indomain] / lam_obs - 1) 
        # y[indomain] = C * Flux * np.exp(-0.5 * exponand * exponand) / (lam_obs * Dispersion * SQRT_2_PI)
        
        return y
        
    # Derivative with respect to every parameter
    def fit_deriv(self, x, Redshift, Flux, Dispersion):
        # 1 + z
        oneplusz = (1 + Redshift)
        # observed center
        lam_obs = oneplusz * self.center
        C2 = C*C
        DDll = lam_obs*lam_obs*Dispersion*Dispersion

        # # Zero outside region
        # indomain = np.logical_and(x > self.domain[0],x < self.domain[1])        
        
        # exponand = ( C / Dispersion ) * ( x[indomain] / lam_obs - 1) 

        # d_Flux = np.zeros(x.shape)    
        # d_Flux[indomain] = np.exp(-0.5 * exponand * exponand) * C / (SQRT_2_PI * Dispersion * lam_obs)
        
        exponand = ( C / Dispersion ) * ( x / lam_obs - 1) 

        d_Flux = np.zeros(x.shape)    
        d_Flux = np.exp(-0.5 * exponand * exponand) * C / (SQRT_2_PI * Dispersion * lam_obs)


        d_Redshift = Flux * d_Flux * \
                     (C2*x*x - C2*x*lam_obs + DDll) / \
                     (oneplusz*DDll)

        d_Dispersion = Flux * d_Flux * \
                        (C2*(lam_obs - x)*(lam_obs - x) - DDll) / \
                        (DDll*Dispersion)
        
        return [d_Redshift, d_Flux, d_Dispersion]

# SSP Continuum
class SSPContinuum(PolynomialModel):

    """
    Continuum from SSPs
    """

    inputs = ('x',)
    outputs = ('y',)
    _separable = False

    def __init__(self, spectrum, n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):

        # Keep track of spectrum
        self.spectrum = spectrum
        self.region = np.invert(self.spectrum.emission_region)

        # List SSPs
        self.ssp_dir = os.path.dirname(os.path.abspath(__file__))+'/SSPs/'
        self.ssp_names = np.sort([x for x in os.listdir(self.ssp_dir) if '.fits' in x])

        # Get SSPs
        self.ssps = []
        for ssp_name in self.ssp_names:
            with fits.open(self.ssp_dir+ssp_name) as f:
                h = f[0].header
                flux = f[0].data
                self.ssps.append(flux)
        self.ssp_wav = (np.arange(flux.size) - h['CRPIX1'] + 1)*h['CDELT1'] + h['CRVAL1']
        self.ssps = np.array(self.ssps)

        # Call Polynomial Class
        super().__init__(
            len(self.ssp_names)+3, n_models=n_models, model_set_axis=model_set_axis,
            name=name, meta=meta, **params)

        # Set bounds
        self.bounds[self.param_names[0]] = (spectrum.z - 0.005,spectrum.z + 0.005)
        self.bounds[self.param_names[1]] = (6/SIGMA_TO_FWHM,20/SIGMA_TO_FWHM)
        for pname in self.param_names[2:]: 
            self.bounds[pname] = (0,None)

        # Set initial parameters
        medians = np.array([np.median(s[s>0]) for s in self.ssps])
        self.parameters = np.append([self.spectrum.z,8.4/SIGMA_TO_FWHM,1.5],np.nanmedian(self.spectrum.flux)/((len(self.ssp_names)+1)*medians))

    def prepare_inputs(self, x, **kwargs):
        inputs, format_info = super().prepare_inputs(x, **kwargs)

        x = inputs[0]
        return (x,), format_info

    def evaluate(self, x, *coeffs):

        coeffs = np.array(coeffs)
        z = coeffs[0]
        sigma = coeffs[1]
        pli = coeffs[2]
        plc = coeffs[3]
        coeffs = coeffs[4:]
        pl = plc*((x/(3000*(1+z))**(pli))
        ssps = self.ssps

        if not self.fixed[self.param_names[0]]:

            # Final number is variance of E-MILES templates
            var_conv = sigma*sigma - 10.616522503600239
            kernel = np.exp(-np.square((self.ssp_wav - 4500))/(2*var_conv))
            ssps = np.array([fftconvolve(ssp,kernel[kernel>0]/kernel.sum(),'same') for ssp in ssps])
            ssps = spectres(self.spectrum.wav,self.ssp_wav*(1+z),ssps)

        return pl + np.dot(coeffs.T,ssps[:,self.region])[0]
        
    def get_names(self):
        return ['Continuum_Redshift','Continuum_Dispersion','Power_Law_Index','Power_Law_Coefficient']+[x.replace('.fits','') for x in self.ssp_names]

    def fix_params(self):

        # Fix params
        for pn in self.param_names[0:2]: self.fixed[pn] = True

        # Final number is variance of E-MILES templates
        var_conv = self.parameters[1]*self.parameters[1] - 10.616522503600239
        kernel = np.exp(-np.square((self.ssp_wav - 4500))/(2*var_conv))
        self.ssps = spectres(self.spectrum.wav,self.ssp_wav*(1+self.parameters[0]),np.array([fftconvolve(ssp,kernel[kernel>0]/kernel.sum(),'same') for ssp in self.ssps]))

    def set_region(self,region):
        self.region = region

# Polynomail Continuum
class PolyContinuum(PolynomialModel):

    r"""
    Copy of 1D Polynomial model with modification.

    It is defined as:

    .. math::

        P = \sum_{i=0}^{i=n}C_{i} * x^{i}

    within the give domain, otherwise 0
    
    Parameters
    ----------
    degree : int
        degree of the series
    domain : list 
    **params : dict
        keyword: value pairs, representing parameter_name: value

    """

    inputs = ('x',)
    outputs = ('y',)
    _separable = True

    def __init__(self, degree, domain, n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
                
        # Set parameters
        self.domain = domain
        self.center = (domain[0] + domain[1])/2
        super().__init__(
            degree, n_models=n_models, model_set_axis=model_set_axis,
            name=name, meta=meta, **params)

    def prepare_inputs(self, x, **kwargs):
        inputs, format_info = super().prepare_inputs(x, **kwargs)

        x = inputs[0]
        return (x,), format_info

    def evaluate(self, x, *coeffs):
        y = self.horner(x - self.center, coeffs)
        y[np.logical_or(x < self.domain[0],x > self.domain[1])] = 0
        return y

    def fit_deriv(self, x, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        params : throw away parameter
            parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        # Zero outside domain
        indomain = np.logical_and(x > self.domain[0],x < self.domain[1])
        v = np.empty((self.degree + 1,) + x.shape, dtype=float)
        v[0] = 0
        v[0][indomain] = 1
        if self.degree > 0:
            v[1] = 0
            v[1][indomain] = x[indomain] - self.center
            for i in range(2, self.degree + 1):
                v[i] = v[i - 1] * x
        return np.rollaxis(v, 0, v.ndim)

    @staticmethod
    def horner(x, coeffs):
        if len(coeffs) == 1:
            c0 = coeffs[-1] * np.ones_like(x, subok=False)
        else:
            c0 = coeffs[-1]
            for i in range(2, len(coeffs) + 1):
                c0 = coeffs[-i] + c0 * x
        return c0