""" Main Function """

# Packages
import numpy as np
import pandas as pd
import multiprocessing as mp
from astropy.table import Table
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

# HQUAILS supporting files
import PARAMS
import RegionFinding as RF
import FittingModel as FM
# emissionLines: Emission Lines Dictionary

# Main Function
def HQUAILS(outname,spectra,emissionLines):

	## Verify Emission Line Dictionary ##
	if not verifyDict(emissionLines):
		print('Unable to verify emission line dictionary, exiting.')
		return
	## Verify Emission Line Dictionary ##

	## Identify fitting regions
	regions_master = RF.identifyComplexes(emissionLines,tol=region_width)
	## Identify fitting regions

	## Prepare Inputs ##
	inputs = []
	names = []
	results = []
	for spectrum in spectra:
		names.append(spectrum[0])
		inputs.append((spectrum,emissionLines,regions_master,background_degree,maxiter,fthresh))
	
		# Not multiprocessing
		results.append(ProcessSpectrum(spectrum,emissionLines,regions_master,\
						background_degree,maxiter,fthresh))
	
	## Prepare Inputs ##
	
	## Multiprocessing ##
	if num_process == None: num_process = mp.cpu_count() - 1
	pool = mp.Pool(processes=num_process)
	results = pool.starmap(ProcessSpectrum, inputs)
	## Multiprocessing ##
	
	## Gather Results ##
	df = pd.concat([result[0] for result in results], axis=0, ignore_index=True, sort=False)

	return df

	df.insert(0,'Object',names)
	t = Table.from_pandas(df)
	t.write(outname,overwrite = True)
	## Gather Results ##

	return t

# Load in spectrum
def LoadSpectrum(filename,z):

	# Import spectrum
	spectrum = pyfits.getdata(filename,1)

	# Only take good values
	weight = spectrum['ivar']
	good = weight > 0
	wav = 10**spectrum['loglam'][good]
	flux = spectrum['flux'][good]
	
	return wav,flux,weight[good],z

# Get fit parameters for galaxy
def ProcessSpectrum(spectrum,emissionLines_master,regions_master,background_degree,maxiter,fthresh):

	# Load Spectrum
	full_spectrum = LoadSpectrum(spectrum[0],spectrum[1])

	## Seting up regions list and emissionLines dictionary for this spectrum ##
	regions,emissionLines =  RF.specRegionAndLines(full_spectrum,emissionLines_master,regions_master) 
	## Seting up regions list and emissionLines dictionary for this spectrum ##

	## Limit spectrum to regions and good values ##
	limited_spectrum = LimitSpectrum(full_spectrum,regions)
	## Limit spectrum to regions and good values ##

	## Model Fitting ##
	model,param_names,cov = FM.FitSpectrum(limited_spectrum,emissionLines,regions,background_degree,maxiter,fthresh)
	## Model Fitting ##

	## Add regions to parameters ##
	for i,region in enumerate(regions):
		param_names.append('Background_'+str(i)+'_Low')
		param_names.append('Background_'+str(i)+'_High')
		parameters = np.concatenate([model.parameters,region])
	## Add regions to parameters ##
	
	## Plotting ##
	figname = spectrum[0].split('/')[-1].split('.')[-2] + '.pdf'
	fig, ax = plt.subplots(figsize = (10,7))
	ax.step(full_spectrum[0],full_spectrum[1],'k')
	ax.step(full_spectrum[0],model(full_spectrum[0]),'r')
	fig.savefig(figname)
	plt.close(fig)
	
	## Plotting ##
	
	return pd.DataFrame([parameters],columns=param_names),cov

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

# Names
names = [np.genfromtxt('SDSSWISEAGN.csv',delimiter=',',dtype='U100,f8',names=['File','z'])[27]]

test = HQUAILS('test.fits',names,emissionLines)
# print(test.columns)