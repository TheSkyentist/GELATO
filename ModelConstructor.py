''' Build Model for Spectrum '''

import numpy as np
import CustomModels as CM
from astropy.modeling import fitting

# Construct Full Model with F-tests for additional parameters
def ConstructModel(spectrum,emissionLines,regions,background_degree,maxiter):

	# Build initial model
	model,param_names = BuildModel(spectrum,emissionLines,regions,background_degree)

	# Tie params
	model = TieParams(model,emissionLines,param_names)
	
	# Fit initial model to the data
	model = FitModel(spectrum,model,emissionLines,maxiter)

	return model,param_names

# Fit Model and remove parameters that aren't fittable
def FitModel(spectrum,init_model,emissionLines,maxiter):
	
	# Unpack
	wav,flux,weight,redshift = spectrum
	
	# Fit model
	fit = fitting.LevMarLSQFitter()
	fit_model = fit(init_model,wav,flux,weights=weight,maxiter=maxiter)
	
	return init_model


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
	for region in regions:
		
		# Name of region
		name = 'Background_' + str((region[0] + region[1])/2) +'_'
		
		# Add model
		model_components.append(CM.ContinuumBackground(background_degree,region))
		
		# Collect param names
		for pname in model_components[-1].param_names:
			param_names.append(name+pname)
	
	# Over all emission lines
	for z in emissionLines.keys():
		for species in emissionLines[z].keys():
			for line in emissionLines[z][species][0]:
				name =  z + '_' + species + '_' + str(line[0]) + '_'
				
				# Add model
				model_components.append(CM.SpectralFeature(center = line[0]))
				
				# Collect param names
				for pname in model_components[-1].param_names:
					param_names.append(name+pname)
	model = np.sum(model_components)
	## Build Base Model

	## Assign Starting Values ##
	for i,region in enumerate(regions):
		# Set intercept to median value
		inregion = np.logical_and(wav > region[0],wav < region[1])
		model.parameters[(background_degree+1)*i] = np.median(flux[inregion])
		
	# Over all emission lines
	for z in emissionLines.keys():
		for species in emissionLines[z].keys():
			for line in emissionLines[z][species][0]:
				# Find parameter name prefix
				name =  z + '_' + species + '_' + str(line[0]) + '_'
				
				# Set redshift to starting redshift
				redshift_index = param_names.index(name + 'Redshift')
				model.parameters[redshift_index] = redshift
				
				# Set height to starting height
				height_index = param_names.index(name + 'Height')
				model.parameters[height_index] = np.max(flux[inregion]) - np.median(flux[inregion])
	## Assign Starting Values ##
	
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
def GenTieFunc(index,scale=1):

	def TieFunc(model):
		return model.parameters[index]*scale

	return TieFunc

