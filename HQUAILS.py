import scipy
import scipy.stats as stats
import numpy as np
from astropy.modeling import models, fitting, custom_model

# Temporary packages
import matplotlib.pyplot as plt
from importlib import reload


# HQUAILS supporting files
import RegionFinding
reload(RegionFinding)
# emissionLines: Emission Lines Dictionary


# Emission line dictionary
emissionLines = {
'Narrow':{
	'OIII':([(5006.77,1),(4958.83,0.350)],1,['Broad']),
	'Hbeta':([(4861.32,1)],1,['Broad']),
	'NII':([(6583.34,1),(6547.96,0.340)],1,['Broad']),
	'Halpha':([(6562.80,1)],1,['Broad'])},
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

	# Verify Emission Line Dictionary
	if not verifyDict(emissionLines):
		print('Unable to verify emission line dictionary, exiting.')
		return
	print('Emission line dictionary verified.')

	## Identify fitting regions
	regions = RegionFinding.identifyComplexes(emissionLines,tol=region_width)
	print('Identififed emission line complexes.')

# Iterate over spectra (eventually)
	wav,flux,weight,z = x,y,weights,0
	

	
	# Check if there is spectral coverage of the regions
	for region in regions:
		# If not...
		if (wav.min() > region[0]*(1+z)) or (wav.max() < region[1]*(1+z)):
			# ...remove region
			regions.remove(region)
			
	## Moodel Constructor
	
	
	print(regions)
# Verify if emission line dictionary is in correct format. 
def verifyDict(emissionLines):

	# Check if Emission Lines are in dictionary
	if not type(emissionLines) == dict:
		print('Emission lines must be in dictionary.')
		return False
	
	# Check if redshift group is in a dictionary
	for z in emissionLines.keys():
		if not type(emissionLines[z]) == dict:
			print('Each redshift must be a dictionary.')
			return False
		if not type(z) == str:
			print('Redshift key must be a string.')
			return False
		
		# Check if species is in an appropriate tuple or list
		for species in emissionLines[z]:
			if not ((type(emissionLines[z][species]) == tuple) or \
					(type(emissionLines[z][species]) == list)):
				print('Each species must be a tuple or list.')
				return False
			if not type(species) == str:
				print('Species key must be a string.')
				return False
			if not len(emissionLines[z][species]) == 3:
				print('Species must be length 3.')
				return False
			# Check if lines are in a list or tuple
			if not ((type(emissionLines[z][species][0]) == list) or \
					(type(emissionLines[z][species][0]) == tuple)):
				print('Lines must be in a list or tuple.')
				return False
				
			# Check if flag is an good int
			if not type(emissionLines[z][species][1]) == int:
				print('Lines flag must be an int.')
				return False
			if not emissionLines[z][species][1] >= 0:
				print('Lines flag must be a positive int.')
				return False
				
			# Check if additional components are in a list 
			if not (type(emissionLines[z][species][2]) == list):
				print('Line component redshift groups must be a list.')
				return False
							
			# Check if each line is an appropriate tuple or list
			for line in emissionLines[z][species][0]:
				if not (type(line) == tuple) or (type(line) == list):
					print('Each line must be a tuple or list.')
					return False
				if not len(line) == 2:
					print('Each line tuple must have two elements.')
					return False
				for l in line:
					if not ((type(l) == float) or (type(l) == int)):
						print('Each element of tuple must be a float or int.')
						return False
			
			# Check if flag bit sum is equal to length of 
			if not (sum([int(bit) for bit in bin(emissionLines[z][species][1])[2:]]) == \
					len(emissionLines[z][species][2])):
				print('Flag number does not match number of additional components.')
				return False
						
	return True

HQUAILS(emissionLines)
	