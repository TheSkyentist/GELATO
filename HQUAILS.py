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

	name = path.split('/')[-1]
	print("HQUAILS running for",name)

	## Load in Spectrum ##
	spectrum = SC.Spectrum(path,z,params)
	print("Loaded spectrum:",name)

	## Create Base Model ##
	model,param_names = BM.BuildModel(spectrum)
	print("Base model created:",name)

	## Fit Additional Components ##
	model,param_names = FM.FitComponents(spectrum,model,param_names)
	print("Additional components added:",name)

	# Bootstrap
	print("Beggining bootstrap (this may take a while):",name)
	parameters = np.array([FM.FitBoot(spectrum,model) for i in range(params['NBoot'])])
	print("Bootstrap iterations finished:",name)

	## Save Results
	pyfits.BinTableHDU.from_columns([pyfits.Column(name=pn,format='D',array=p) for p,pn in zip(parameters.T,param_names)]).writeto(params['OutFolder']+name.replace('.fits','-results.fits'),overwrite=True)
	print("Results saved:",name)

	## Plotting ##
	if params['Plotting']:
		# Set model parameters to median values
		model.parameters = np.median(parameters,0)
		PL.Plot(spectrum,model,path)
		print("Figure saved:",name)

	print("HQUAILS finished running on":,name)

def header():
	print("Welcome to HQUAILS")
	print("Handy QUAsar emIssion Line fitS (pronounced like 'quails')")
	print("Developed by R. E. Hviding")
	