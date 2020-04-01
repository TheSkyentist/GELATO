""" Spectrum Class """

# Packages
import copy
import numpy as np
from astropy.io import fits

# Constants
C = 299792.458 # km/s

class Spectrum:
    
    def __init__(self,path,z,p):

        # Object Parameter List
        self.p = p

        # Object's redshift
        self.z = z

        # Load Spectrum
        spectrum = fits.getdata(path,1)

        # Only take good values
        weight = spectrum['ivar']
        good = weight > 0

        # Initial data
        self.wav = 10**spectrum['loglam'][good]
        self.flux = spectrum['flux'][good]
        self.weight = weight[good]

        # Create regions and lines
        self.regionAndLines()
        self.reduceRegions()
        self.LimitSpectrum()
    
    # Return reduced regions and emission lines based on spectrum wavelength
    def regionAndLines(self):

        # Find regions and emission line lists
        self.regions = []
        for group in self.p['EmissionGroups']:
            for species in group['Species']:
                for line in species['Lines']:
                    linewav = line['Wavelength']*(1+self.z)
                    
                    # Ensure there is spectral coverage of the line
                    linewidth = linewav*self.p['LineDataWidth']/(2*C)
                    if np.logical_and(np.any(self.wav > linewav + linewidth),np.any(self.wav < linewav - linewidth)):
                        dellam = linewav*self.p['RegionWidth']/(2*C)
                        self.regions.append([np.max([linewav - dellam,self.wav[0]]), np.min([linewav + dellam,self.wav[-1]])])

        # Remove emission lines not in regions
        # For each emission line
        eG = copy.deepcopy(self.p['EmissionGroups'])
        for group in eG:
            gname = group['Name']
            for species in group['Species']:
                sname = species['Name']
                for line in species['Lines']:
                    linewav = line['Wavelength']

                    # Check for removal
                    remove = True # Initialize as removing
                    for region in self.regions: # Check if it is in any region
                        if (region[0] < linewav*(1+self.z)) and (linewav*(1+self.z) < region[1]):
                            remove = False # If it is, set as remove 
                    
                    # Remove if necessary
                    if remove:
                        for g in self.p['EmissionGroups']:
                            for s in g['Species']:
                                for l in s['Lines']:
                                    # Check if the same then remove
                                    if ((g['Name'] == gname) and (s['Name'] == sname) and (l['Wavelength'] == linewav)):
                                        s['Lines'].remove(l)
                                # If species is empty, delete it
                                if (len(s['Lines']) == 0):
                                    g['Species'].remove(s)

    # Reduce number of regions so none overlap
    def reduceRegions(self):

        # Sort regions
        self.regions = np.sort(self.regions,0)

        # Initalize loop
        for i,region in enumerate(self.regions[:-1]):
            if self.regions[i+1][0] < region[1]:
                self.regions[i] = [region[0],self.regions[i+1][1]]
                self.regions = np.delete(self.regions,i+1,0)
                self.reduceRegions()
                break
    
    # Pair down spectrum
    def LimitSpectrum(self):

        # Only take data in regions
        inregion = []
        for region in self.regions:
            inregion.append(np.logical_and(region[0]<self.wav,self.wav<region[1]))
        inregion = np.logical_or.reduce(inregion)
        self.wav = self.wav[inregion]
        self.flux = self.flux[inregion]
        self.weight = self.weight[inregion]
        self.sqrtweight = np.sqrt(self.weight)
        self.sigma = 1/self.sqrtweight

    # Boostrap the Flux
    def Boostrap(self):

        return np.random.normal(self.flux,self.sigma)