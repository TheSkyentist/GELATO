import scipy
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting, custom_model

# Plotting
plt.close('all')

# Random seed
np.random.seed(42)

# Chi Squared
def chi2(model,x,y,weights):
	r = y - model(x)
	return np.sum(r*r*weights)

# Two Models
g1 = models.Gaussian1D(amplitude=2, mean=0, stddev=0.2)
g2 = models.Gaussian1D(amplitude=2.5, mean=0, stddev=0.05)

# x data
N 		= 200
x 		= np.linspace(-1, 1, N)
sigmas 	= np.random.uniform(0.4,0.5, x.shape)
weights = 1/np.square(sigmas)
y_true 	= g1(x) + g2(x)
y 		= g1(x) + g2(x) + np.random.normal(0., sigmas, x.shape)
# y 		= g1(x) + np.random.normal(0., sigmas, x.shape)

# Fittings
fit 	= fitting.LevMarLSQFitter()
g1 = models.Gaussian1D(amplitude=2, mean=0, stddev=0.1)
g2 = models.Gaussian1D(amplitude=2.5, mean=0, stddev=0.1)
g_redu_fit	= fit(model = g1, x = x, y = y, weights = weights)
g_full_fit	= fit(model = g1+g2, x = x, y = y, weights = weights)

# Number of parameters
df_redu = N - g_redu_fit.parameters.size
df_full = N - g_full_fit.parameters.size

# Calculate Chi2
RSS_redu = chi2(g_redu_fit,x,y,weights)
RSS_full = chi2(g_full_fit,x,y,weights)

# F-test
F_distrib 	= scipy.stats.f(df_redu - df_full,df_full)
# F_distrib 	= scipy.stats.f(df_full,df_redu - df_full)
F_value		= ((RSS_redu - RSS_full)/(df_redu - df_full))/(RSS_full/df_full)
plt.plot(np.linspace(0,5,200),F_distrib.pdf(np.linspace(0.001,5,200)))

plt.show()

print(F_distrib.cdf(F_value))

# Plotting
plt.errorbar(x,y,yerr=sigmas)
plt.errorbar(x,y_true)

plt.plot(x,g1(x)+g2(x),label='Guess')
plt.plot(x,g_redu_fit(x),label = 'Reduced')
plt.plot(x,g_full_fit(x),label = 'Full')
plt.legend()

plt.show()