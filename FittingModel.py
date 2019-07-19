''' Build Model for Spectrum '''

import copy
import numpy as np
from astropy.modeling import fitting

import CustomModels as CM
import AdditionalComponents as AC
import FischerTest as FT

# Construct Full Model with F-tests for additional parameters
def FitSpectrum(spectrum,emissionLines_base,regions,background_degree,maxiter,threshold):
	
	# Build initial model
	base_model,base_param_names = BuildModel(spectrum,emissionLines_base,regions,background_degree)
	# Fit model
	base_model,_ = FitModel(spectrum,base_model,maxiter)

	# Find number of flags
	flags = 0
	for group in emissionLines_base.keys():
		for species in emissionLines_base[group].keys():
			flagbits = bin(emissionLines_base[group][species][1])[2:]
			flags += len(flagbits.replace('0',''))
					
	# Initial F-tests
	accepted = []
	accepted_models = []
	accepted_names = []
	for i in range(flags):
	
		# Add new component
		emissionLines = AddComplexity(emissionLines_base,i)
		model,param_names = BuildModel(spectrum,emissionLines,regions,background_degree)
		
		# Source starting parameters
		model = SourceParams(model,param_names,base_model,base_param_names)		
		
		# Fit model
		model,_ = FitModel(spectrum,model,maxiter)

		if FT.FTest(spectrum,regions,base_model,model,threshold):
			accepted.append(i)
			accepted_models.append(model)
			accepted_names.append(param_names)

	# Add new component
	emissionLines = AddComplexity(emissionLines_base,accepted)
	model,param_names = BuildModel(spectrum,emissionLines,regions,background_degree)
	
	# Source starting parameters
	for source,source_params in zip(accepted_models,accepted_names):
		model = SourceParams(model,param_names,source,source_params)
	model = SourceParams(model,param_names,base_model,base_param_names)		
	
	# Fit model
	model,cov = FitModel(spectrum,model,maxiter)		
		
	return model.parameters,param_names,cov

# Fit Model
def FitModel(spectrum,init_model,maxiter):
	
	# Unpack
	wav,flux,weight,redshift = spectrum
	
	# Fit model
	fit = fitting.LevMarLSQFitter()
	fit_model = fit(init_model,wav,flux,weights=weight,maxiter=maxiter)
	
	return fit_model,fit.fit_info['param_cov']

# Carry over parameters from another model
def SourceParams(model,param_names,source,source_params):
	
	# Iterate over source_oaram names
	for i,source_param in enumerate(source_params):
		if source_param in param_names:
			index = param_names.index(source_param)
			model.parameters[index] = source.parameters[i]
			
	return model

# Build the model from the emissionLines and spectrum with initial guess
def BuildModel(spectrum,emissionLines,regions,background_degree):

	# Unpack
	wav,flux,_,redshift = spectrum 
	
	## Build Base Model
	# True param names
	param_names	= []
	
	# Model
	model_components = []
	
	# Add background regions
	for i,region in enumerate(regions):
		
		# Name of region
		name = 'Background_' + str(i) +'_'
		
		# Generate model
		background = CM.ContinuumBackground(background_degree,region)
		
		# Add starting parameters
		inregion = np.logical_and(wav > region[0],wav < region[1])
		background.parameters[0] = np.median(flux[inregion])
		
		# Add model
		model_components.append(background)
		
		# Collect param names
		for pname in model_components[-1].param_names:
			param_names.append(name+pname)
	
	# Over all emission lines
	for z in emissionLines.keys():
		for species in emissionLines[z].keys():				
		
			for line in emissionLines[z][species][0]:
				name =  z + '_' + species + '_' + str(line[0]) + '_'
				
				# If additional component
				if emissionLines[z][species][1] < 0:
					model = AC.AddComponent(emissionLines[z][species][1],line[0],spectrum, regions)
					
				# IF original model component
				else:
					model = CM.SpectralFeature(center = line[0],spectrum = spectrum,regions = regions)

				# Add model
				model_components.append(model)
				
				# Collect param names
				for pname in model_components[-1].param_names:
					param_names.append(name+pname)
	model = np.sum(model_components)
	## Build Base Model

	## Tie parameters ##
	model = TieParams(model,emissionLines,param_names)
	## Tie parameters ##

	return model,param_names

# Tied all model parameters
def	TieParams(model,emissionLines,param_names):

	## Tie parameters ##
	for z in emissionLines.keys():
		reference_z = True
		for species in emissionLines[z].keys():
				
			reference_line = True
			for line in emissionLines[z][species][0]:

				# Find parameter name prefix
				name =  z + '_' + species + '_' + str(line[0]) + '_'
				
				## Tie Redshift
				# Check for first line in redshift and make tie function
				if reference_z:
					reference_z = False
					TieRedshift = GenTieFunc(param_names.index(name+'Redshift'))
				# Otherwise tie redshift
				else:
					redshift_index = param_names.index(name+'Redshift')
					model.tied[model.param_names[redshift_index]] = TieRedshift
					
				## Tie Sigma and Flux/Height
				# Check for first line in species and make tie function or initialize scale
				if reference_line:
					# Sigma
					reference_line	= False
					TieSigma 		= GenTieFunc(param_names.index(name+'Sigma'))
					
					# Flux
					reference_flux 	= line[1]
					index_flux		= param_names.index(name+'Height')

				# Otherwise tie param
				else:
					# Sigma
					sigma_index = param_names.index(name+'Sigma')
					model.tied[model.param_names[sigma_index]] = TieSigma				
					
					# Redshhift
					height_index = param_names.index(name+'Height')
					TieHeight = GenTieFunc(index_flux,scale=line[1]/reference_flux)
					model.tied[model.param_names[height_index]] = TieHeight				
				## Tie Tie Sigma and Flux/Height
	##Tie parameters ##

	return model
	
# Generate Tie Paramater Function, needed to preserve static values in function
# Note: tried using lambda function in place, but it didn't work. 
def GenTieFunc(index,scale=1):

	def TieFunc(model):
		return model.parameters[index]*scale

	return TieFunc

# Add additional component to a model
def AddComplexity(emissionLines,index):
	
	# If multiple indices, add all of them
	if hasattr(index,'__iter__'):
		for i in index:
			emissionLines = AddComplexity(emissionLines,i)
		return emissionLines

	# Deepcopy
	emissionLines = copy.deepcopy(emissionLines)
			
	# Keeping track
	i = 0
	
	# Iterate in emission line dictionary 
	for z in emissionLines.keys():
		for species in emissionLines[z].keys():
			
			# Flag
			flag 		= bin(emissionLines[z][species][1])[2:]
			flag_len 	= np.sum([int(bit) for bit in flag])
			
			# If we have a flag
			if (flag_len > 0):
				
				# Check that our index is in the range
				if (index >= i) and (index < i + flag_len):		
				
					# Position in flagged bits
					j = 0
						
					# Iterate over bits in flag
					for bit in flag:
						
						# If we arrive at the correct index
						if bit == '1':
						
							# If we've arrived at the index
							if index == i:
							
								# Construct the added entry
								entry = (emissionLines[z][species][0],int('-0b1'+j*'0',2),[])
							
								# Redshift key
								z = emissionLines[z][species][2][j]
								# Add it if it doesn't exist
								if z not in emissionLines.keys():
									emissionLines[z] = {}
								
								# Species key
								species = species + '-' + AC.ComponentName(j)

								emissionLines[z][species] = entry
								
								return emissionLines
						
							# Increment along the bit							
							i += 1
							j += 1
					
				i += flag_len