""" Main Function """

# Packages
import copy
import argparse
import numpy as np
# import pandas as pd
import multiprocessing as mp
from functools import partial
# from astropy.table import Table
# import astropy.io.fits as pyfits

# HQUAILS supporting files
# import Plotting as PL
# import FittingModel as FM
# import RegionFinding as RF
import SpectrumClass as SC

# Get fit parameters for galaxy
# outfolder,spectrum,emissionLines_master,regions_master,background_degree,maxiter,fthresh
def ProcessSpectrum(obj,p):

	### Load in Spectrum ###
	spectrum = SC.Spectrum(obj,p)

	return
	## Seting up regions list and emissionLines dictionary for this spectrum ##
	regions,emissionLines =  RF.specRegionAndLines(full_spectrum,emissionLines_master,regions_master) 
	## Seting up regions list and emissionLines dictionary for this spectrum ##

	## Limit spectrum to regions and good values ##
	limited_spectrum = LimitSpectrum(full_spectrum,regions)
	## Limit spectrum to regions and good values ##

	## Model Fitting ##
	model,param_names,parameters = FM.FitSpectrum(limited_spectrum,emissionLines,regions,background_degree,maxiter,fthresh,n_boot)
	sigma = np.std(parameters,0)
	parameters = np.median(parameters,0)
	model.parameters = parameters
	## Model Fitting ##

	## Create ouput ##
	out = np.empty(parameters.size*2)
	out[0::2] = parameters
	out[1::2] = sigma
	outnames = []
	for p in param_names:
		outnames.append(p)
		outnames.append(p+'_sig')
	## Create ouput ##

	## Add regions to parameters ##
	for i,region in enumerate(regions):
		outnames = np.append(outnames,'Background_'+str(i)+'_Low')
		outnames = np.append(outnames,'Background_'+str(i)+'_High')
		out = np.concatenate([out,region])
	## Add regions to parameters ##

	## Plotting ##
	PL.Plot(outfolder,spectrum[0],model,full_spectrum,regions)
	## Plotting ##

	print('Finished fitting '+spectrum[0].split('/')[-1])
	return pd.DataFrame(data=out.reshape((1,out.size)),columns=outnames)

# Pair down spectrum
def LimitSpectrum(spectrum,regions):

	# Unpack
	wav,flux,weight,z = spectrum 
	
	# Only take data in regions
	inregion = []
	for region in regions:
		inregion.append(np.logical_and(region[0]<wav,wav<region[1]))
	inregion = np.logical_or.reduce(inregion)
	wav = wav[inregion]
	flux = flux[inregion]
	weight = weight[inregion]
	
	return 	wav,flux,weight,z

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

# Main Function
if __name__ == "__main__":

	## Parse Arguements to find Parameter File ##
	parser = argparse.ArgumentParser()
	parser.add_argument("Settings", type=str, help="Path to settings file")
	args = parser.parse_args()
	exec(open(args.Settings).read())
	p = PARAMETERS()
	## Parse Arguements to find Parameter File ##

	## Verify Emission Line Dictionary ##
	if not verifyDict(p['emissionLines']):
		print('Unable to verify emission line dictionary, exiting.')
		sys.exit(1)
	## Verify Emission Line Dictionary ##

	# ## Prepare Inputs ##
	object_list = p.pop('names')
	if p['num_process'] == 1: # Single Thread
		results = [ProcessSpectrum(obj,copy.deepcopy(p)) for obj in object_list]
	else: # Multi-threading
		if p['num_process'] == None: p['num_process'] = mp.cpu_count() - 1
		pool = mp.Pool(processes=p['num_process'])
		results = pool.map(partial(ProcessSpectrum, p = p),object_list)
		
	# ## Process Spectra ##
	
	# ## Gather and Write Results ##
	# df = pd.concat([result for result in results], ignore_index=True, sort=False)
	# df.insert(0,'Object',names)
	# t = Table.from_pandas(df)
	# t.write(p['outfolder']+'results.fits',overwrite = True)
	# ## Gather and Write Results ##