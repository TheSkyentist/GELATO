''' Additional components '''

import numpy as np
import CustomModels as CM

SQRT_2_PI = np.sqrt(2*np.pi)

def ComponentName(index):

	'''
	For each bit position, what is the name given to each additional component such that
	they can be identified in the list of model name parameters.
	Can be as many as desired as long as paired with a component below. 
	'''

	if index == 0:
		return 'Broad'

	if index == 1:
		return 'Absorption'

def AddComponent(flag, line, spectrum, regions):

	# Unpack
	wav,flux,weight,redshift = spectrum 

	'''
	For each bit position, what is the model that we implement,
	and how do we set the initial parameters.
	Must return astropy model with edits
	'''
	
	# Find index
	flag 	= bin(flag)[3:]
	index 	= len(flag) - flag.index('1') - 1
	
	
	''' Broad '''
	if index == 0: 
	
		# Broad line model
		# Use wider default dispersion
		model = CM.SpectralFeature(center = line,spectrum = spectrum,regions = regions, Dispersion = 1200)
			
		return model

	''' Absorption '''
	if index == 1: 
	
		# Absorption line model
		# Use wider default dispersion		
		model = CM.SpectralFeature(center = line,spectrum = spectrum,regions = regions, Dispersion = 600)
		model.Flux.bounds = (None,0) # Must be non-positive
		
		# Set default depth
		for region in regions:
			if (line*(1+redshift) < region[1]) and (line*(1+redshift) > region[0]):
				break
		inregion = np.logical_and(wav > region[0],wav < region[1])				
		model.Flux = ( np.min(flux[inregion]) - np.median(flux[inregion]) ) * model.Dispersion * SQRT_2_PI
				
		return model