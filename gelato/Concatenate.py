""" Concatenate Results """

# Import packages
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack

# Concatenate results
def concatfromresults(p,objects):

    if p["Verbose"]:
        print("Combining gelato...")

    # Initalize loop
    first = True
    N = 10000 # Number of objects we will concatenate at once
    i = 0 # Index
    while i*N < len(objects):
        
        tables = [] # Initilize tables list
        paths = objects['File'][i*N:(i+1)*N] # Get subsample

        # Iterate over results
        for path in paths:
            
            # Load name and parameters
            name = path.split('/')[-1].replace('.fits','')
            parameters = fits.getdata(p['OutFolder']+name+'-results.fits')
            
            # Initalize Lists
            data = [name]
            names = ['Name']
            dtype = [np.unicode_] + [np.float_ for i in range(2*len(parameters.columns.names))]
            
            # Iterate over columns and add
            for n in parameters.columns.names:
                
                ps = parameters[n]
                ps = ps[np.invert(np.isinf(ps))]

                # Add medians
                data.append(np.nanmedian(ps))
                names.append(n)

                # Add errors
                data.append(np.nanstd(ps))
                names.append(n+'_err')

            tables.append(Table(data = np.array(data), names = names,dtype=dtype))

        table = vstack(tables,join_type = 'outer')

        # If first entry, save to disk directly
        if first:
            table.write(p['OutFolder']+'GELATO-results.fits',overwrite=True)
            first = False

        # Otherise load table and append to it
        else:
            results = Table.read(p['OutFolder']+'GELATO-results.fits')
            vstack([results,table],join_type = 'outer').write(p['OutFolder']+'GELATO-results.fits',overwrite=True)

        i += 1

    if p["Verbose"]:
        print("Gelato combined.")