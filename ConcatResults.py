#! /usr/bin/env python

""" Concatenate Results """

# Import packages
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack

# Concatenate results
def concatfromresults(p,objects):

    if p["Verbose"]:
        print("Concatenating Results...")

    # Initalize list of tables
    first = True
    for path in objects['File']:

        # Load name and parameters
        name = path.split('/')[-1].replace('.fits','')
        parameters = fits.getdata(p['OutFolder']+name+'-results.fits')

        # Initalize Lists
        data = [name]
        names = ['Name']
        dtype = [np.unicode_] + [np.float_ for i in range(2*len(parameters.columns.names))]

        # Iterate over columns and add
        for n in parameters.columns.names:
            
            # Add medians
            data.append(np.median(parameters[n]))
            names.append(n)

            # Add errors
            data.append(np.std(parameters[n]))
            names.append(n+'_err')

        # Make table
        table = Table(data = np.array(data), names = names,dtype=dtype)

        # If first entry, initialized table
        if first:
            table.write(p['OutFolder']+'HQUAILS-results.fits',overwrite=True)
            first = False
        # Otherise load table and append to it
        else:
            results = Table.read(p['OutFolder']+'HQUAILS-results.fits')
            vstack([results,table],join_type = 'outer').write(p['OutFolder']+'HQUAILS-results.fits',overwrite=True)

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