
import numpy as np
import pandas as pd
import multiprocessing as mp
from astropy.modeling import fitting
import matplotlib.pyplot as plt

def FitResults(spectrum,model,maxiter,param_names,num_boostrap,num_process):

	# Assign number of pools
	if num_process == None:
		num_process = mp.cpu_count() - 1
		
	# Generate input tuples
	inputs = []
	for i in range(num_boostrap):
		inputs.append((spectrum,model,maxiter))

	results = []
	for input in inputs:
		results.append(BootstrapFit(input))

	# Multiprocessing
# 	pool = mp.Pool(processes=num_process)
# 	results = pool.map(BootstrapFit, inputs)
# 	pool.close()
# 	pool.join()
# 	

	results = np.array(results)[:,4]
	plt.hist(results)
	plt.show()
	
	

# Fit Model
def BootstrapFit(input):

	spectrum,model,maxiter = input
	
	# Unpack
	wav,flux,weight,redshift = spectrum
	
	# Resample
	flux = np.random.normal(loc=flux,scale=1/np.sqrt(weight))
	
	# Fit model
	fit = fitting.LevMarLSQFitter()
	fit_model = fit(model,wav,flux,weights=weight,maxiter=maxiter)
	
	print(fit.fit_info)
	# Put output
	return fit_model.parameters