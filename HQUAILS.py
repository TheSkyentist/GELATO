""" Main Function """

# Temporary packages
import matplotlib.pyplot as plt
plt.close('all')
from astropy.modeling import models, fitting, custom_model

import numpy as np
import pandas as pd
import multiprocessing as mp
from astropy.table import Table

# HQUAILS supporting files
import RegionFinding as RF
import FittingModel as FM
# emissionLines: Emission Lines Dictionary


# Emission line dictionary
# Full testing for chip gaps
emissionLines = {
'AGN':{
	'OIII':([(5006.77,1),(4958.83,0.350)],1,['Broad']),	
	'NII':([(6583.34,1),(6547.96,0.340)],1,['Broad'])},
'Galaxy':{
	'Halpha':([(6562.80,1)],1,['Broad']),
	'Hbeta':([(4861.32,1)],3,['Broad','Absorption']),
	'Htest':([(6500,1)],0,[])
		},
'Broad':{},
'Absorption':{},
'Bad':{
	'Test':([(6000,1)],1,['Broad'])
	},
}

emissionLines = {
'AGN':{
	'OIII':([(5006.77,1),(4958.83,0.350)],1,['Outflow']),	
	'NII':([(6583.34,1),(6547.96,0.340)],1,['Broad'])},
'Galaxy':{
	'Halpha':([(6562.80,1)],1,['Broad']),
	'Hbeta':([(4861.32,1)],3,['Broad','Galaxy'])
		},
'Bad':{},
}

emissionLines = {
'AGN':{
	'OIII':([(5006.77,1),(4958.83,0.350)],1,['Outflow']),	
	'NII':([(6583.34,1),(6547.96,0.340)],0,[])},
'Galaxy':{
	'Halpha':([(6562.80,1)],1,['Broad']),
	'Hbeta':([(4861.32,1)],3,['Broad','Galaxy'])
		},
'Bad':{},
}

# Fake data
# Random seed
np.random.seed(100)
		
## Fake spectrum
z = 0.5
true_x = np.arange(7000,10100,1)
x = true_x#np.concatenate([np.arange(7000,8950,1),np.arange(9050,9715,1),np.arange(9790,10000,1)]) # Angstroms

# Continuum
y = 0.5*np.ones(x.shape) # Flux
y += (x - 8250.0)*.0001 + -(x-8250.0)*(x-8250.0)*.0000001
# Balmer
y += models.Gaussian1D(amplitude=1, mean=(1+z)*6562.80, stddev=5)(x)
y += models.Gaussian1D(amplitude=2/2.86, mean=(1+z)*4861.32, stddev=5)(x)
y += models.Gaussian1D(amplitude=1, mean=(1+z)*6562.80, stddev=20)(x)
y += models.Gaussian1D(amplitude=-0.4, mean=(1+z)*4861.32, stddev=10)(x)
# NII
y += models.Gaussian1D(amplitude=1, mean=(1+z)*6583.34, stddev=5)(x)
y += models.Gaussian1D(amplitude=0.34, mean=(1+z)*6547.96, stddev=5)(x)
# OIII
y += models.Gaussian1D(amplitude=0.5, mean=(1+z)*5006.77, stddev=5)(x)
y += models.Gaussian1D(amplitude=0.35/2, mean=(1+z)*4958.83, stddev=5)(x)
y += models.Gaussian1D(amplitude=0.5, mean=(1+z-0.002)*5006.77, stddev=20)(x)
y += models.Gaussian1D(amplitude=0.35/2, mean=(1+z-0.002)*4958.83, stddev=20)(x)

# Test lines
# y += models.Gaussian1D(amplitude=0.5, mean=(1+z)*6500, stddev=5)(x)
# y += models.Gaussian1D(amplitude=0.5, mean=(1+z)*6000, stddev=5)(x)

# Add error
sigmas = np.random.uniform(0.01,0.05, x.shape)
weights = 1/np.square(sigmas)
y += np.random.normal(0., sigmas, x.shape)

spectrum1 = ['test1',x,y,weights,z]
spectrum2 = ['test2',x,y,weights,z]

# Main Function
def HQUAILS(outname,spectra,emissionLines,region_width=100,background_degree=1,maxiter=1000,fthresh=0.95,num_process=None):

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
	for spectrum in spectra:
		names.append(spectrum[0])
		inputs.append((spectrum[1:],emissionLines,regions_master,background_degree,maxiter,fthresh))
	## Prepare Inputs ##
	
	## Multiprocessing ##
	if num_process == None: num_process = mp.cpu_count() -1
	pool = mp.Pool(processes=num_process)
	results = pool.starmap(ProcessSpectrum, inputs)
	## Multiprocessing ##
	
	## Gather Results ##
	df = pd.concat([result[0] for result in results], axis=0, ignore_index=True)
	df.insert(0,'Object',names)
	t = Table.from_pandas(df)
	t.write(outname,overwrite = True)
	## Gather Results ##

	return t

def ProcessSpectrum(spectrum,emissionLines_master,regions_master,background_degree,maxiter,fthresh):

	## Seting up regions list and emissionLines dictionary for this spectrum ##
	regions,emissionLines =  RF.specRegionAndLines(spectrum,emissionLines_master,regions_master) 
	## Seting up regions list and emissionLines dictionary for this spectrum ##

	## Limit spectrum to regions and good values ##
	spectrum = LimitSpectrum(spectrum,regions)
	## Limit spectrum to regions and good values ##

	## Model Fitting ##
	parameters,param_names,cov = FM.FitSpectrum(spectrum,emissionLines,regions,background_degree,maxiter,fthresh)
	## Model Fitting ##

	return pd.DataFrame([parameters],columns=param_names),cov

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

def LimitSpectrum(spectrum,regions):

	# Unpack
	wav,flux,weight,z = spectrum 
	
	# Only take good values
	good = weight > 0
	wav = wav[good]
	flux = flux[good]
	weight = weight[good]
	
	# Only take data in regions
	inregion = []
	for region in regions:
		inregion.append(np.logical_and(region[0]<wav,wav<region[1]))
	inregion = np.logical_or.reduce(inregion)
	wav = wav[inregion]
	flux = flux[inregion]
	weight = weight[inregion]
	
	return 	wav,flux,weight,z

print(HQUAILS('test.fits',[spectrum1,spectrum2], emissionLines))