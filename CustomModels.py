import numpy as np
from astropy.modeling import Fittable1DModel,Parameter

from astropy.modeling import fitting
import matplotlib.pyplot as plt
plt.close('all')

# Global vars
FLOAT_EPSILON = float(np.finfo(np.float32).tiny)
GAUSSIAN_SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))
SQRT_2_PI = np.sqrt(2*np.pi)

## Custom Astropy Modeling Models
# Spectral Feature Class (Gaussian)
class SpectralFeature(Fittable1DModel):

	# Parameters
	redshift = Parameter(default=0)
	height	 = Parameter(default=0)
	sigma	 = Parameter(default=1,bounds=(FLOAT_EPSILON, None)) # Must be positive

	def __init__(self, center = 0, redshift = redshift.default, height = height.default,
				sigma = sigma.default, **kwargs):
		self.center = center
		super().__init__(
			redshift = redshift, height=height, sigma=sigma, **kwargs)
		
	@property
	def fwhm(self):
		"""Gaussian full sigma at half maximum."""
		return self.sigma * GAUSSIAN_SIGMA_TO_FWHM

	@property
	def flux(self):
		"""Gaussian flux."""
		return self.height * self.sigma * SQRT_2_PI

	def evaluate(self, x, redshift, height, sigma):
		"""
		Gaussian1D model function.
		"""
		exponand = (x - self.center*(1 + redshift)) / sigma
		return height * np.exp(- 0.5 * exponand * exponand)
		
	def fit_deriv(self, x, redshift, height, sigma):
		
		exponand = (x - self.center*(1 + redshift)) / sigma
		
		d_height 	= np.exp(-0.5 * exponand * exponand)
		d_redshift 	= height * d_height * exponand / sigma
		d_sigma 	= height * d_height * exponand * exponand / sigma
		
		return [d_redshift, d_height, d_sigma]

# Spectral Feature Class (Gaussian)
class ContinuumBackground(Fittable1DModel):

	# Default _param_names list; this will be filled in by the __init__
	_param_names = ()

	@property
	def param_names(self):
		"""Coefficient names generated based on the model's polynomial degree
		and number of dimensions.
			
		Subclasses should implement this to return parameter names in the
		desired format.
		
		On most `Model` classes this is a class attribute, but for polynomial
		models it is an instance attribute since each polynomial model instance
		can have different parameters depending on the degree of the polynomial
		and the number of dimensions, for example.
		"""

		return self._param_names
        
	def __getattr__(self, attr):
		if self._param_names and attr in self._param_names:
			return Parameter(attr, default=0.0, model=self)
	
        raise AttributeError(attr)
        
    def __setattr__(self, attr, value):
        # TODO: Support a means of specifying default values for coefficients
        # Check for self._ndim first--if it hasn't been defined then the
        # instance hasn't been initialized yet and self.param_names probably
        # won't work.
        # This has to vaguely duplicate the functionality of
        # Parameter.__set__.
        # TODO: I wonder if there might be a way around that though...
        if attr[0] != '_' and self._param_names and attr in self._param_names:
            param = Parameter(attr, default=0.0, model=self)
            # This is a little hackish, but we can actually reuse the
            # Parameter.__set__ method here
            param.__set__(self, value)
        else:
            super().__setattr__(attr, value)
        
np.random.seed(42)
x = np.arange(-50,50,0.1)
g_true = SpectralFeature(center = 1, redshift = 10.1, height = 10, sigma = 5)
y = g_true(x) + np.random.normal(loc=0,scale=1,size=x.shape)

g_init = SpectralFeature(center = 1, redshift = 10, height = 1, sigma = 10)

fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, x, y)

plt.plot(x,y)
plt.plot(x,g(x))
plt.show()


