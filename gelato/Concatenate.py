""" Concatenate Results """

# Import packages
import copy
import numpy as np
from os import path
from tqdm import tqdm
from astropy.io import fits
from astropy.table import Table,vstack

# Order names correctly
def orderNames(table):

    # Order of prefixes
    order = ['Name','SSP_Redshift','SSP_','PowerLaw']

    # Get columns
    colnames = copy.deepcopy(table.colnames)

    # Make empty columns
    names = []

    # Iterate over prefixes
    for o in order:
        newnames = np.sort([n for n in colnames if o in n]).tolist()
        for n in newnames: colnames.remove(n)
        names += newnames

    # Get final names
    lastnames = np.sort([n for n in colnames if 'rChi2' in n]).tolist()
    for n in lastnames: colnames.remove(n) 

    # Order everything else
    names += np.sort(colnames).tolist() + lastnames

    return names

# Concatenate results
def concatfromresults(p,objects):

    if p["Verbose"]:
        print("Combining gelato...")

    # Initalize loop
    first = True
    N = 10000 # Number of objects we will concatenate at once
    i = 0 # Index
    pbar = tqdm(total=len(objects),disable=not p["Verbose"]) # Progress Bar
    
    while i*N < len(objects):

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
            
            # Append to list
            tables.append(Table(data=np.array(data),names=names,dtype=dtype))

        # Stack tables, include all columns
        table = vstack(tables,join_type = 'outer')
        if not type(table.mask) == type(None): # NaN for missing values
            for c in table.colnames: table[c][table.mask[c]] = np.nan

        # If first entry, save to disk directly
        out = path.join(p['OutFolder'],'GELATO-results.fits')
        if first:
            table[orderNames(table)].write(out,overwrite=True)
            first = False

        # Otherise load table and append to it
        else:
            results = vstack([Table.read(out),table],join_type='outer')
            results[orderNames(results)].write(out,overwrite=True)

        # Update Progress Bar
        pbar.update(len(spaths))

        i += 1

    # Close progress bar
    pbar.close()

    if p["Verbose"]:
        print("Gelato combined.")