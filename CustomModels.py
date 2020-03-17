""" Custom Astropy Models """

from astropy.modeling import Fittable1DModel,Parameter
from astropy.modeling.polynomial import PolynomialModel
import numpy as np

# Global vars
FLOAT_EPSILON = float(np.finfo(np.float32).tiny)
GAUSSIAN_Sigma_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))
SQRT_2_PI = np.sqrt(2*np.pi)
C = 299792.458 # km/s

## Custom Astropy Modeling Models
# Spectral Feature Class (Gaussian)
class SpectralFeature(Fittable1DModel):

    """
    Very similar to the 1D Gaussian but parametrized differently. 
    """

    # Parameters
    Redshift     = Parameter(default=0)
    Flux         = Parameter(default=0) # Must be non-negative for emission
    Dispersion    = Parameter(default=150,bounds=(75,500)) # Reasonable range

    def __init__(self, center, spectrum, Dispersion = Dispersion.default, **kwargs):
        
        # Find starting height
        # Find corresponding region
        for region in spectrum.regions:
            if (center*(1+spectrum.z) < region[1]) and (center*(1+spectrum.z) > region[0]):
                break
        inregion = np.logical_and(spectrum.wav > region[0],spectrum.wav < region[1])
        Height = np.abs(spectrum.flux[np.argmin(np.abs(spectrum.wav - center*(1+spectrum.z)))] - np.median(spectrum.flux[inregion]))
        Flux = Height * (1 + spectrum.z) * Dispersion * center * SQRT_2_PI / C
        
        # Set parameters
        self.center = center
        self.domain = region
        super().__init__(Redshift = spectrum.z, Flux=Flux, Dispersion=Dispersion, **kwargs)
        self.Redshift.bounds = (spectrum.z - 0.01,spectrum.z + 0.01)
        self.Flux.bounds = (0,1.5*Flux*self.Dispersion.bounds[1]/self.Dispersion) # Set positive bounds

    @property
    def sigma(self):
        """Gaussian full Sigma at half maximum."""
        return self.Dispersion * self.center / C 
        
    @property
    def fwhm(self):
        """Gaussian full Sigma at half maximum."""
        return self.sigma() * GAUSSIAN_Sigma_TO_FWHM

    def evaluate(self, x, Redshift, Flux, Dispersion):
        """
        Gaussian1D model function.
        """
        
        # (1 + z) * lambda
        lam_obs = (1 + Redshift) * self.center
        
        # Zero outside region
        indomain = np.logical_and(x > self.domain[0],x < self.domain[1])        
        
        # Assign values
        y = np.zeros(x.shape)
        exponand = (C / Dispersion) * (x[indomain] / lam_obs - 1) 
        y[indomain] = C * Flux * np.exp(-0.5 * exponand * exponand) / (lam_obs * Dispersion * SQRT_2_PI)
        
        return y
        
    # Derivative with respect to every parameter
    def fit_deriv(self, x, Redshift, Flux, Dispersion):
        # 1 + z
        oneplusz = (1 + Redshift)
        # observed center
        lam_obs = oneplusz * self.center
        C2 = C*C
        DDll = lam_obs*lam_obs*Dispersion*Dispersion

        # Zero outside region
        indomain = np.logical_and(x > self.domain[0],x < self.domain[1])        
        
        exponand = ( C / Dispersion ) * ( x[indomain] / lam_obs - 1) 

        d_Flux = np.zeros(x.shape)    
        d_Flux[indomain] = np.exp(-0.5 * exponand * exponand) * C / (SQRT_2_PI * Dispersion * lam_obs)

        d_Redshift = Flux * d_Flux * \
                     (C2*x*x - C2*x*lam_obs + DDll) / \
                     (oneplusz*DDll)

        d_Dispersion = Flux * d_Flux * \
                        (C2*(lam_obs - x)*(lam_obs - x) - DDll) / \
                        (DDll*Dispersion)
        
        return [d_Redshift, d_Flux, d_Dispersion]

# Spectral Feature Class (Gaussian)
class ContinuumBackground(PolynomialModel):

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

    n_inputs = 1
    n_outputs = 1
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

    @property
    def input_units(self):
        if self.degree == 0 or self.c1.unit is None:
            return None
        else:
            return {'x': self.c0.unit / self.c1.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        mapping = []
        for i in range(self.degree + 1):
            par = getattr(self, 'c{0}'.format(i))
            mapping.append((par.name, outputs_unit['y'] / inputs_unit['x'] ** i))
        return OrderedDict(mapping)