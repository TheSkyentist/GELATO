""" Statistical F-Test Package """
import numpy as np
import scipy.stats as stats

# Preform an f-test
def FTest(model1,x1,model2,x2,spectrum,args):

    # Find full and reduced model
    if (model1.nparams() < model2.nparams()):
        model_full = model2
        x_full = x2
        model_redu = model1
        x_redu = x1
    else: 
        model_full = model1
        x_full = x1
        model_redu = model2
        x_full = x2

    # Degrees of freedom
    N = args[0].size

    df_redu = N - model_redu.nparams()
    df_full = N - model_full.nparams()

    # Calculate Chi2
    RSS_redu = Chi2(model_redu,x_redu,args)
    RSS_full = Chi2(model_full,x_full,args)

    # F-test
    F_value = ((RSS_redu - RSS_full)/(df_redu - df_full))/(RSS_full/df_full)
    F_distrib = stats.f(df_redu - df_full,df_full)

    # Greater than threshold?
    return F_distrib.cdf(F_value) > spectrum.p['FThresh']

# Calculate Akaike Information Criterion
def AIC(model,p,args):

    # Correct up to a constant
    return Chi2(model,p,args) + 2*model.nparams()

def Chi2(model,p,args):

    return np.square(model.residual(p,*args)).sum()