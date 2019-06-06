""" Statistical F-Test Package """

import numpy as np
import scipy.stats as stats

# Temporary packages
import matplotlib.pyplot as plt
plt.close('all')

# Preform an f-test
def FTest(spectrum,regions,model1,model2,threshold):

	# Unpack
	wav,flux,weight,redshift = spectrum

# 	plt.step(wav,flux,where='mid')
# 	plt.step(wav,model1(wav),where='mid')
# 	plt.step(wav,model2(wav),where='mid')
# 	plt.show()
	
	# Find full and reduced model
	if (model1.parameters.size < model2.parameters.size):
		model_full = model2
		model_redu = model1
	else: 
		model_full = model1
		model_redu = model2
	
	# Degrees of freedom
	N = 0
	for region in regions:
		N += np.sum(np.logical_and(region[0] < wav,wav < region[1]))
	df_redu = N - model_redu.parameters.size
	df_full = N - model_full.parameters.size	

	# Calculate Chi2
	RSS_redu = Chi2(model_redu,wav,flux,weight)
	RSS_full = Chi2(model_full,wav,flux,weight)

	# F-test
	F_value		= ((RSS_redu - RSS_full)/(df_redu - df_full))/(RSS_full/df_full)
	F_distrib 	= stats.f(df_redu - df_full,df_full)

    # Greater than threshold?
	return F_distrib.cdf(F_value) > threshold

# Chi Squared of model
def Chi2(model,x,y,weights):
	r = y - model(x) # Residual
	return np.sum(r*r*weights)