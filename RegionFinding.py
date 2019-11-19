""" Region Finding Code """

import copy
import numpy as np 

# Reduce number of regions so none overlap
def reduceRegions(regions):

	# Sort regions
	regions = np.sort(regions,0)

	# Initalize loop
	for i,region in enumerate(regions[:-1]):
		if regions[i+1][0] < region[1]:
			regions[i] = [region[0],regions[i+1][1]]
			regions = reduceRegions(np.delete(regions,i+1,0))
			break

	# Return regions			
	return regions

# Identify line complexes present with given emission lines	
def identifyComplexes(emissionLines,tol=50):

	# Initialize list of line centers
	regions = []

	# Iterate over redshifts
	for group in emissionLines.keys():
		# Iterate over species
		for species in emissionLines[group].keys():
			# Iterate over line centers
			for line in emissionLines[group][species][0]:
		
				# Append Line Centers
				regions.append([line[0]-tol,line[0]+tol])
			
	return reduceRegions(regions)
	
# Return reduced regions and emission lines based on spectrum wavelength
def specRegionAndLines(spectrum,emissionLines_master,regions_master):
	# Unpack
	wav,_,_,z = spectrum 
	
	# Deep copy
	regions 		= copy.deepcopy(regions_master)
	regions			= [[region[0]*(1+z),region[1]*(1+z)] for region in regions]
	emissionLines 	= copy.deepcopy(emissionLines_master)
	
	# Check if there is spectral coverage of the regions
	for region in regions:
		# If not...
		if np.sum(np.logical_and(wav > region[0],wav < region[1])) == 0:
			# ...remove region
			regions.remove(region)
		
	# Remove emission lines not in regions
	# For each emission line	
	for group in emissionLines_master.keys():
		for species in emissionLines_master[group].keys():
			for line in emissionLines_master[group][species][0]:

				# Check for removal
				remove = True # Initialize as removing
				for region in regions: # Check if it is in any region
					if (region[0] < line[0]*(1+z)) and (line[0]*(1+z) < region[1]):
						remove = False # If it is, set as remove 
				
				# Remove if necessary
				if remove:
					emissionLines[group][species][0].remove(line)
			
			# If species is empty, delete it
			if (len(emissionLines[group][species][0]) == 0):
				del emissionLines[group][species]		

	return regions,emissionLines