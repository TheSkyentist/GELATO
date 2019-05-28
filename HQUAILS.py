""" Main Function """

# Temporary packages
import matplotlib.pyplot as plt
plt.close('all')
from importlib import reload
import numpy as np
from astropy.modeling import models, fitting, custom_model

# HQUAILS supporting files
import RegionFinding as RF
import ModelConstructor as MC
# emissionLines: Emission Lines Dictionary


# Emission line dictionary
emissionLines = {
'AGN':{
	'OIII':([(5006.77,1),(4958.83,0.350)],1,['Broad']),
	
	'NII':([(6583.34,1),(6547.96,0.340)],1,['Broad'])},
'Galaxy':{
	'Halpha':([(6562.80,1)],1,['Broad']),
	'Hbeta':([(4861.32,1)],3,['Broad','Absorption'])
		},
 'Broad':{},
 'Absorption':{}
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
def background(x,slope=0,intercept=1,rightedge=0,leftedge=0):
	center 	= (rightedge + leftedge)/2
	y 		= np.zeros(x.shape)
	y[np.logical_and(x > rightedge,x < leftedge)] = slope*(x - center) + intercept
	
## Fake spectrum
z = 0.5
x = np.arange(6000,10500,1) # Angstroms
y = 0.5*np.ones(x.shape) # Flux
# Balmer
# y += models.Gaussian1D(amplitude=0.5, mean=(1+z)*6562.80, stddev=20)(x)
y += models.Gaussian1D(amplitude=0.5, mean=(1+z)*6562.80, stddev=5)(x)
y += models.Gaussian1D(amplitude=1/2.86, mean=(1+z)*4861.32, stddev=5)(x)
# y += models.Gaussian1D(amplitude=-0.1, mean=(1+z)*4861.32, stddev=10)(x)
# NII
y += models.Gaussian1D(amplitude=1, mean=(1+z)*6583.34, stddev=5)(x)
y += models.Gaussian1D(amplitude=0.34, mean=(1+z)*6547.96, stddev=5)(x)
# OIII
# y += models.Gaussian1D(amplitude=0.5, mean=(1+z-0.002)*5006.77, stddev=20)(x)
# y += models.Gaussian1D(amplitude=0.35, mean=(1+z-0.002)*4958.83, stddev=20)(x)
y += models.Gaussian1D(amplitude=0.5, mean=(1+z)*5006.77, stddev=5)(x)
y += models.Gaussian1D(amplitude=0.35, mean=(1+z)*4958.83, stddev=5)(x)
# Add error
sigmas = np.random.uniform(0.01,0.05, x.shape)
weights = 1/np.square(sigmas)
y += np.random.normal(0., sigmas, x.shape)
y += (x - 8250.0)*.0001 + -(x-8250.0)*(x-8250.0)*.0000001

spectrum = [x,y,weights,z]

# Main Function
def HQUAILS(emissionLines_master,region_width=50,background_degree=1):

	# Verify Emission Line Dictionary
	if not verifyDict(emissionLines_master):
		print('Unable to verify emission line dictionary, exiting.')
		return
	print('Emission line dictionary verified.')

	## Identify fitting regions
	regions_master = RF.identifyComplexes(emissionLines_master,tol=region_width)
	print('Identififed emission line complexes.')

	''' Iterate over spectra (eventually) '''
	wav,flux,weight,z = spectrum 
	
	## Seting up regions list and emissionLines dictionary for this spectrum ##
	regions,emissionLines =  RF.specRegionAndLines(spectrum,emissionLines_master,regions_master) 
	## Seting up regions list and emissionLines dictionary for this spectrum ##

	## Model construction ##
	model,param_names = MC.ConstructModel(spectrum,emissionLines,regions,background_degree)
	
	
	return model
	## Model construction ##
	
			

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
			if not (type(species) == str):
				print('Species key must be a string.')
				return False
			if not (len(emissionLines[z][species]) == 3):
				print('Species must be length 3.')
				return False
			# Check if lines are in a list or tuple
			if not (type(emissionLines[z][species][0]) == list):
				print('Lines must be in a list.')
				return False
				
			# Check if flag is an good int
			if not (type(emissionLines[z][species][1]) == int):
				print('Lines flag must be an int.')
				return False
			if not (emissionLines[z][species][1] >= 0):
				print('Lines flag must be a positive int.')
				return False
				
			# Check if additional components are in a list 
			if not (type(emissionLines[z][species][2]) == list):
				print('Line component redshift groups must be a list.')
				return False
			for component in emissionLines[z][species][2]:
				if not component in emissionLines.keys():
					print('Additional component redshift must be in emission dictionary.')
					return False
					
			# Check if each line is an appropriate tuple or list
			for line in emissionLines[z][species][0]:
				if not (type(line) == tuple) or (type(line) == list):
					print('Each line must be a tuple or list.')
					return False
				if not (len(line) == 2):
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

wav,flux,weight,z = spectrum 

test = HQUAILS(emissionLines)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(test,wav,flux,weights=weight)

plt.plot(wav,flux)
# plt.plot(wav,g(wav))
plt.show()
