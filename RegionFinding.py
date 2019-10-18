""" Region Finding Code """

import copy
import numpy as np

# If two regions overlap, return true
def overlap(a,b):
	
	return min(a[1], b[1]) - max(a[0], b[0]) > 0

# Reduce number of regions so none overlap
def reduceRegions(regions):
	
	# Check every unique region pairing 
	for i in range(len(regions)):
	
		# Break loop if region list shorted
		if (i >= len(regions)):
			break

		for j in range(i+1,len(regions)):
		
			# Break loop if region list shorted
			if (i >= len(regions)) or (j >= len(regions)):
				break
			
			# If they overlap
			if overlap(regions[i],regions[j]):
							
				# Append compound region to list
				regions.append([min(regions[i][0],regions[j][0]),max(regions[i][1],regions[j][1])])
				
				# Remove parent regions
				regions.remove(regions[j]) # Remove j first as j > i
				regions.remove(regions[i])
				
				# Recurse
				reduceRegions(regions)
	
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