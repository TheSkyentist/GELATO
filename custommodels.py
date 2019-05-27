import scipy
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting, custom_model

# Plotting
plt.close('all')

# Random seed
np.random.seed(42)

# Chi Squared
def chi2(model,x,y,weights):
	r = y - model(x)
	return np.sum(r*r*weights)
	
# Define emission line custom model
@custom_model
def emissionLine(x,redshift=0,height=0,sigma=1,center=0):
	exponand = (x - center*(1 + redshift)) # What goes in the exponent
	return height * np.exp( -0.5 * exponand * exponand)

# Define linear background custom model
@custom_model
def background(x,slope,intercept,rightedge,leftedge):
	center 	= (rightedge + leftedge)/2
	y 		= np.zeros(x.shape)
	y[np.logical_and(x > rightedge,x < leftedge)] = slope*(x - center) + intercept
	
## Fake spectrum
x = np.arange(4000,7000,1) # Angstroms
y = 0.5*np.ones(x.shape) # Flux
# Balmer
y += models.Gaussian1D(amplitude=1, mean=6562.80, stddev=5)(x)
y += models.Gaussian1D(amplitude=1/2.86, mean=4861.32, stddev=5)(x)
# NII
y += models.Gaussian1D(amplitude=1, mean=6583.34, stddev=5)(x)
y += models.Gaussian1D(amplitude=0.34, mean=6547.96, stddev=5)(x)
# OIII
y += models.Gaussian1D(amplitude=1, mean=5006.77, stddev=5)(x)
y += models.Gaussian1D(amplitude=0.35, mean=4958.83, stddev=5)(x)
# Add error
sigmas = np.random.uniform(0.01,0.05, x.shape)
weights = 1/np.square(sigmas)
y += np.random.normal(0., sigmas, x.shape)

emissionLines = {
'Narrow':{
	'OIII':([(5006.77,1),(4958.83,0.350)],0),
	'Hbeta':([(4861.32,1)],0),
	'NII':([(6583.34,1),(6547.96,0.340)],0),
	'Halpha':([(6562.80,1)],0)},
'Broad':{
	'OIII':([(5006.77,1),(4958.83,0.350)],0),
	'Hbeta':([(4861.32,1)],0),
	'NII':([(6583.34,1),(6547.96,0.340)],0),
	'Halpha':([(6562.80,1)],0)}
}

plt.step(x,y,where='mid')
plt.show()