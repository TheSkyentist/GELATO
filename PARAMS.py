'''
    Parameter file 
'''

# Import packages
import numpy as np

# Names and redshifts
names = np.genfromtxt('../WISEAGN/SDSSWISEAGN.csv',delimiter=',',dtype='U100,f8',names=['File','z'])[27:30]

# Output name
outfolder = '../WISEAGN/Results/' # Folder where Results will go

# Emission line dictionary
emissionLines = {
    'AGN':{
        'OIII':([(5006.77,1),(4958.83,0.350)],1,['Outflow']),	
        'NII':([(6583.34,1),(6547.96,0.340)],0,[])},
    'Galaxy':{
        'Halpha':([(6562.80,1)],1,['Broad']),
        'Hbeta':([(4861.32,1)],3,['Broad','Galaxy'])
            },
    }

# Fitting Parameters
region_width        = 100   # Width of region around emission lines (AA)
background_degree   = 1     # Degree of background 
maxiter             = 1000  # Maximum LM iterations
n_boot              = 500   # Number of Bootstrap Iterations
fthresh             = 0.95  # Probability threshold for f-test
num_process         = 2     # Number of threads, None uses all threads except 1