""" Main Function """

# Packages
import numpy as np
import astropy.io.fits as pyfits

# HQUAILS supporting files
import Plotting as PL
import BuildModel as BM
import FittingModel as FM
import SpectrumClass as SC

# Get fit parameters for spectrum
def HQUAILS(params,path,z):

	## Load in Spectrum ##
	spectrum = SC.Spectrum(path,z,params)

	## Create Base Model ##
	model,param_names = BM.BuildModel(spectrum)

	## Fit Additional Components ##
	model,parameters,param_names = FM.FitComponents(spectrum,model,param_names)

	## Add in regions ##
	for i,region in enumerate(spectrum.regions):
		param_names = np.append(param_names,'Background_'+str(i)+'_Low')
		param_names = np.append(param_names,'Background_'+str(i)+'_High')
		parameters = np.hstack((parameters,region[0]*np.ones((parameters.shape[0],1))))
		parameters = np.hstack((parameters,region[0]*np.ones((parameters.shape[0],1))))

	## Save Results
	pyfits.BinTableHDU.from_columns([pyfits.Column(name=pn,format='D',array=p) for p,pn in zip(parameters.T,param_names)]).writeto(params['OutFolder']+path.split('/')[-1].replace('.fits','-results.fits'),overwrite=True)

	## Plotting ##
	if params['Plotting']:
		PL.Plot(spectrum,model,path)

