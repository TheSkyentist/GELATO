''' Additional components '''

import numpy as np
import CustomModels as CM

SQRT_2_PI = np.sqrt(2*np.pi)

def ComponentName(index):

	'''
	For each bit position, from right to left, what is the name given to each additional component such that they can be identified in the list of model name parameters.
	So the rightmost bit corresponds to index 0.
	Can be as many as desired as long as paired with a component below. 
	'''

	if index == 0:
		return 'Broad'

	if index == 1:
		return 'Absorption'

def AddComponent(flag, line, spectrum):

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
		model = CM.SpectralFeature(center = line,spectrum = spectrum, Dispersion = 400)
			
		return model

	''' Absorption '''
	if index == 1: 
	
		# Absorption line model
		# Use wider default dispersion		
		model = CM.SpectralFeature(center = line,spectrum = spectrum, Dispersion = 600)
		model.Flux.bounds = (None,0) # Must be non-positive
		
		# Set default depth
		for region in spectrum.regions:
			if (line*(1+spectrum.z) < region[1]) and (line*(1+spectrum.z) > region[0]):
				break
		inregion = np.logical_and(spectrum.wav > region[0],spectrum.wav < region[1])				
		model.Flux = ( np.min(spectrum.flux[inregion]) - np.median(spectrum.flux[inregion]) ) * model.Dispersion * SQRT_2_PI
				
		return model