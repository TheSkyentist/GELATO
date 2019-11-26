""" Concatenate Results """

# Import packages
import numpy as np
import pandas as pd
import astropy.io.fits as pyfits
from astropy.table import Table,vstack

# Concatenate results
def concatfromresults(p,objects):
    
    tables = []
    for path in objects['File']:
        parameters = pyfits.getdata(p['OutFolder']+path.split('/')[-1].replace('.fits','-results.fits'),1)
        tables.append(Table(data = np.array([np.median(parameters[n]) for n in parameters.columns.names]), names = parameters.names))
    vstack(tables,join_type = 'outer').write(p['OutFolder']+'results.fits')


# Main Function
if __name__ == "__main__":

    import argparse
    import numpy as np
    import ConstructParams as CP

    ## Parse Arguements ##
    parser = argparse.ArgumentParser()
    parser.add_argument('Parameters', type=str, help='Path to parameters file')
    parser.add_argument('ObjectList', type=str, help='Path to object list with paths to spectra and their redshifts.')
    args = parser.parse_args()
    p = CP.construct(args.Parameters)
    ## Parse Arguements

    ## Assemble Objects
    objects = np.genfromtxt(args.ObjectList,delimiter=',',dtype='U100,f8',names=['File','z'])
    ## Assemble Objects

    ## Concatenate Results
    concatfromresults(p,objects)