import scipy
import scipy.stats as stats
import numpy as np
from astropy.modeling import models, fitting, custom_model

import matplotlib.pyplot as plt


# HQUAILS supporting files
from RegionFinding import identifyComplexes

# emissionLines: Emission Lines Dictionary


# Emission line dictionary
emissionLines = {
'Narrow':{
	'OIII':([(5006.77,1),(4958.83,0.350)],1,('Broad')),
	'Hbeta':([(4861.32,1)],1,('Broad')),
	'NII':([(6583.34,1),(6547.96,0.340)],1,('Broad')),
	'Halpha':([(6562.80,1)],1,('Broad'))},
'Broad':{}
}

# Fake data
# Random seed
np.random.seed(42)
	
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



# Main Function
def HQUAILS(emissionLines,region_width=50):

	## Iterate over spectra (eventually)
	wav,flux,weight,z = x,y,weights,0
	
	## Identify fitting regions
	regions = identifyComplexes(emissionLines,tol=region_width)
	
	# Check if there is spectral coverage of the regions
	for region in regions:
		# If not...
		if (wav.min() > region[0]*(1+z)) or (wav.max() < region[1]*(1+z)):
			# ...remove region
			regions.remove(region)
			
	## Moodel Constructor
	
	
	print(regions)
		

HQUAILS(emissionLines)
	