""" Custom Astropy Models """

from astropy.modeling import Fittable1DModel,Parameter
from astropy.modeling.polynomial import PolynomialModel
import numpy as np

# Global vars
FLOAT_EPSILON = float(np.finfo(np.float32).tiny)
GAUSSIAN_Sigma_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))
SQRT_2_PI = np.sqrt(2*np.pi)

## Custom Astropy Modeling Models
# Spectral Feature Class (Gaussian)
class SpectralFeature(Fittable1DModel):

	r"""
	Very similar to the 1D Gaussian but parametrized differently. 
	"""

	# Parameters
	Redshift = Parameter(default=0)
	Height	 = Parameter(default=0,bounds=(0,None)) # Must be non-negative for emission
	Sigma	 = Parameter(default=5,bounds=(FLOAT_EPSILON, None)) # Must be positive

	def __init__(self, center = 0, Redshift = Redshift.default, Height = Height.default,
				Sigma = Sigma.default, **kwargs):
		self.center = center
		super().__init__(
			Redshift = Redshift, Height=Height, Sigma=Sigma, **kwargs)
		
	@property
	def fwhm(self):
		"""Gaussian full Sigma at half maximum."""
		return self.Sigma * GAUSSIAN_Sigma_TO_FWHM

	@property
	def flux(self):
		"""Gaussian flux."""
		return self.Height * self.Sigma * SQRT_2_PI

	def evaluate(self, x, Redshift, Height, Sigma):
		"""
		Gaussian1D model function.
		"""
		exponand = (x - self.center*(1 + Redshift)) / Sigma
		return Height * np.exp(- 0.5 * exponand * exponand)
		
	# Derivative with respect to every parameter
	def fit_deriv(self, x, Redshift, Height, Sigma):
		
		exponand = (x - self.center*(1 + Redshift)) / Sigma
		
		d_Height 	= np.exp(-0.5 * exponand * exponand)
		d_Redshift 	= Height * d_Height * exponand / Sigma
		d_Sigma 	= Height * d_Height * exponand * exponand / Sigma
		
		return [d_Redshift, d_Height, d_Sigma]

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

	inputs = ('x',)
	outputs = ('y',)
	_separable = True

	def __init__(self, degree, domain, n_models=None,
				 model_set_axis=None, name=None, meta=None, **params):
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