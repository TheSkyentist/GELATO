"""Spectrum Class"""

import numpy as np

import astropy.io.fits as pyfits

class Spectrum:
    
    def __init__(self,obj,p):

        # Object Parameter List
        self.params = p

        # Object's redshift
        self.redshift = obj[1]

        # Load Spectrum
        spectrum = pyfits.getdata(obj[0],1)

        # Only take good values
        weight = spectrum['ivar']
        good = weight > 0

        # Initial data
        self.wav = 10**spectrum['loglam'][good]
        self.flux = spectrum['flux'][good]
        self.weight = weight[good]

        # Find regions and emission line lists
        self.regions = []
        for group in p['emissionLines'].keys():
            for species in p['emissionLines'][group].keys():
                for line in p['emissionLines'][group][species][0]:
                    self.regions.append([line[0]-p['region_width'],line[0]+p['region_width']])
        self.specRegionAndLines()
        
        
    
    # Return reduced regions and emission lines based on spectrum wavelength
    def specRegionAndLines(self):

        # Check if there is spectral coverage of the regions
        for region in self.regions:
            # If not...
            if np.sum(np.logical_and(self.wav > region[0],self.wav < region[1])) == 0:
                # ...remove region
                self.regions.remove(region)
            
        # # Remove emission lines not in regions
        # # For each emission line	
        # for group in emissionLines_master.keys():
        #     for species in emissionLines_master[group].keys():
        #         for line in emissionLines_master[group][species][0]:

        #             # Check for removal
        #             remove = True # Initialize as removing
        #             for region in regions: # Check if it is in any region
        #                 if (region[0] < line[0]*(1+z)) and (line[0]*(1+z) < region[1]):
        #                     remove = False # If it is, set as remove 
                    
        #             # Remove if necessary
        #             if remove:
        #                 emissionLines[group][species][0].remove(line)
                
        #         # If species is empty, delete it
        #         if (len(emissionLines[group][species][0]) == 0):
        #             del emissionLines[group][species]		

        # return regions,emissionLines