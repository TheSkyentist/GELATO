
import numpy as np
import pandas as pd
import multiprocessing as mp
from astropy.modeling import fitting

# 
print(mp.cpu_count())


def FitResults(spectrum,model,param_names,num_boostrap,num_process):

	# Assign number of pools
	if num_process == None:
		num_process = mp.cpu_count() - 1
		
	pool = mp.Pool(processes=num_pools)
	results = pool.map(BootstrapFit, range(1,7))
	pool.close()
	pool.join()

# Fit Model
def BootstrapFit(spectrum,model,maxiter, output):
	
	# Unpack
	wav,flux,weight,redshift = spectrum
	
	# Resample
	flux = np.random.normal(loc=flux,scale=1/np.sqrt(weight))
	
	# Fit model
	fit = fitting.LevMarLSQFitter()
	fit_model = fit(model,wav,flux,weights=weight,maxiter=maxiter)
	
	# Put output
	return fit_model.parameters