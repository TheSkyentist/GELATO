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
	for z in emissionLines.keys():
		# Iterate over species
		for species in emissionLines[z].keys():
			# Iterate over line centers
			for line in emissionLines[z][species][0]:
		
				# Append Line Centers
				regions.append([line[0]-tol,line[0]+tol])
			
	return reduceRegions(regions)