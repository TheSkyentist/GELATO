""" Statistical F-Test Package """
import numpy as np
import scipy.stats as stats

# Preform an f-test
def FTest(spectrum,model1,model2,fitregion = None):

    if type(fitregion) == type(None): 
        fitregion = np.ones(spectrum.wav.shape,dtype=bool)

    # Find full and reduced model
    if (model1.parameters.size < model2.parameters.size):
        model_full = model2
        model_redu = model1
    else: 
        model_full = model1
        model_redu = model2

    # Degrees of freedom
    # N = 0
    # for region in spectrum.regions:
    #     N += np.sum(np.logical_and(region[0] < spectrum.wav[fitregion],spectrum.wav[fitregion] < region[1]))
    N = fitregion.sum()

    df_redu = N - np.sum([val == False for val in model_redu.tied.values()])
    df_full = N - np.sum([val == False for val in model_full.tied.values()])

    # Calculate Chi2
    RSS_redu = Chi2(model_redu,spectrum.wav[fitregion],spectrum.flux[fitregion],spectrum.weight[fitregion])
    RSS_full = Chi2(model_full,spectrum.wav[fitregion],spectrum.flux[fitregion],spectrum.weight[fitregion])

    # F-test
    F_value = ((RSS_redu - RSS_full)/(df_redu - df_full))/(RSS_full/df_full)
    F_distrib = stats.f(df_redu - df_full,df_full)

    # Greater than threshold?
    return F_distrib.cdf(F_value) > spectrum.p['FThresh']

# Calculate Akaike Information Criterion
def AIC(model,spectrum):

    # Model parameters
    k = np.sum([val == False for val in model.tied.values()])

    # Correct up to a constant
    return Chi2(model,spectrum.wav,spectrum.flux,spectrum.weight) + 2*k

# Chi Squared of model
def Chi2(model,x,y,weights):
    r = y - model(x) # Residual
    return np.sum(r*r*weights)

# Reduce Chi Squared of model
def rChi2(spectrum,model):
    # Degrees of freedom
    N = np.sum([np.sum(np.logical_and(region[0] < spectrum.wav,spectrum.wav < region[1])) for region in spectrum.regions])
    dof = N - np.sum([val == False for val in model.tied.values()])
    return Chi2(model,spectrum.wav,spectrum.flux,spectrum.weight)/dof