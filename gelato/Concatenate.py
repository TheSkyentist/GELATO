""" Concatenate Results """

# Import packages
import numpy as np
from os import path
from astropy.io import fits
from astropy.table import Table,vstack

# Gelato dependecies
import gelato.Utility as U

# Concatenate results
def concatfromresults(p,objects):

    if p["Verbose"]:
        print("Combining gelato...")

    # Initalize loop
    first = True
    N = 10000 # Number of objects we will concatenate at once
    i = 0 # Index
    while i*N < len(objects):

        # Loading bar
        Nbar = 40
        if p['Verbose']: 
                pc = int(100*i*N/len(objects)) # Percentage
                l = int(Nbar*i*N/len(objects)) # Length of bar
                if pc == 0:
                    print('Progress: |'+Nbar*'-'+'|   0%',end='\r')
                elif pc < 10:
                    print('Progress: |'+l*'#'+(Nbar-l)*'-'+'|   '+str(pc)+'%',end='\r')
                else:
                    print('Progress: |'+l*'#'+(Nbar-l)*'-'+'|  '+str(pc)+'%',end='\r')

        tables = [] # Initilize tables list
        spaths = objects['Path'][i*N:(i+1)*N] # Get subsample

        # Iterate over results
        for spath in spaths:
            
            # Load name, output file, and parameters
            name = path.split(spath)[-1]
            if (name[-5:] == '.fits'):
                spath = name[:-5]+'-results.fits'
            elif (name[-8:] == '.fits.gz'):
                spath = name[:-8]+'-results.fits'
            else:
                spath = name+'-results.fits'
            spath = path.join(p['OutFolder'],spath)
            if not path.exists(spath):     
                continue # If doesn't exist, continue
            parameters = fits.getdata(spath,'PARAMS')

            # Initalize Lists
            data = [name]
            names = ['Name']
            dtype = [np.unicode_] + [np.float64 for i in range(2*len(parameters.columns.names))]
            
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
        if not type(table.mask) == type(None):
            for c in table.colnames: table[c][table.mask[c]] = np.nan

        # If first entry, save to disk directly
        out = path.join(p['OutFolder'],'GELATO-results.fits')
        if first:
            table.write(out,overwrite=True)
            first = False

        # Otherise load table and append to it
        else:
            results = Table.read(out)
            vstack([results,table],join_type = 'outer').write(out,overwrite=True)

        i += 1

    if p['Verbose']: print('Progress: |'+Nbar*'#'+'| 100%')

    if p["Verbose"]:
        print("Gelato combined.")