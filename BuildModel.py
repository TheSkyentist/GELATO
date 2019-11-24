"""Build Model for Spectrum"""

# Packages
import numpy as np

# HQUAILS supporting files
import CustomModels as CM
import AdditionalComponents as AC

# Build the model from the emissionLines and spectrum with initial guess
def BuildModel(spectrum, emissionLines=None):
	
	## Build Base Model
	# True param names
	param_names	= []
	
	# Model
	model_components = []
		
	# Add background regions
	for i,region in enumerate(spectrum.regions):
		
		# Name of region
		name = 'Background_' + str(i) +'_'
		
		# Generate model
		background = CM.ContinuumBackground(spectrum.p['BackgroundDeg'],region)
		
		# Add starting parameters
		inregion = np.logical_and(spectrum.wav > region[0],spectrum.wav < region[1])
		background.parameters[0] = np.median(spectrum.flux[inregion])
		
		# Add model
		model_components.append(background)
		
		# Collect param names
		for pname in model_components[-1].param_names:
			param_names.append(name+pname)

	# Check if we were passed an emission lines
	if emissionLines == None:
		emissionLines = spectrum.p['EmissionLines']

	# Over all emission lines
	for z in emissionLines.keys():
		for species in emissionLines[z].keys():				
		
			for line in emissionLines[z][species][0]:
				name =  z + '_' + species + '_' + str(line[0]) + '_'
				
				# If additional component
				if emissionLines[z][species][1] < 0:
					model = AC.AddComponent(emissionLines[z][species][1],line[0],spectrum)
					
				# IF original model component
				else:
					model = CM.SpectralFeature(center = line[0],spectrum = spectrum)

				# Add model
				model_components.append(model)
				
				# Collect param names
				for pname in model_components[-1].param_names:
					param_names.append(name+pname)
	model = np.sum(model_components)
	## Build Base Model

	## Tie parameters ##
	model = TieParams(spectrum,model,param_names,emissionLines)
	## Tie parameters ##

	return model, param_names

# Tied all model parameters
def	TieParams(spectrum,model,param_names,emissionLines):

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
					
				## Tie Dispersion and Flux
				# Check for first line in species and make tie function or initialize scale
				if reference_line:
					# Dispersion
					reference_line	= False
					TieDispersion 		= GenTieFunc(param_names.index(name+'Dispersion'))
					
					# Flux
					reference_flux 	= line[1]
					index_flux		= param_names.index(name+'Flux')

				# Otherwise tie param
				else:
					# Dispersion
					dispersion_index = param_names.index(name+'Dispersion')
					model.tied[model.param_names[dispersion_index]] = TieDispersion				
					
					# Redshhift
					height_index = param_names.index(name+'Flux')
					TieFlux = GenTieFunc(index_flux,scale=line[1]/reference_flux)
					model.tied[model.param_names[height_index]] = TieFlux				
				## Tie Tie Dispersion and Flux
	##Tie parameters ##

	return model

# Generate Tie Paramater Function, needed to preserve static values in function
# Note: tried using lambda function in place, but it didn't work. 
def GenTieFunc(index,scale=1):

	def TieFunc(model):
		return model.parameters[index]*scale

	return TieFunc