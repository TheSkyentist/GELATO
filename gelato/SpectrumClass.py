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
        spectrum = fits.getdata(path)

        # Only take good values
        weight = spectrum['ivar']
        good = weight > 0

        # Initial data
        self.wav = 10**spectrum['loglam'][good]
        self.flux = spectrum['flux'][good]
        self.weight = weight[good]
        self.sqrtweight = np.sqrt(self.weight)
        self.sigma = 1/self.sqrtweight

        # Create regions and lines
        self.regionAndLines()
        self.reduceRegions()
        self.emission_region = np.zeros(self.wav.shape,dtype=bool)
        for r in self.regions:
            self.emission_region[np.logical_and(r[0]<self.wav,self.wav<r[1])] = True
        
        # self.LimitSpectrum()
    
    # Return reduced regions and emission lines based on spectrum wavelength
    def regionAndLines(self):

        # Find regions and emission line lists
        self.regions = []
        eG = copy.deepcopy(self.p['EmissionGroups'])
        for group in eG:
            for species in group['Species']:
                for line in species['Lines']:
                    linewav = line['Wavelength']*(1+self.z)

                    # Ensure there is spectral coverage of the line
                    linewidth = linewav*self.p['LineRegion']/(2*C)
                    if (np.any(self.wav > linewav + linewidth) and np.any(self.wav < linewav - linewidth) and (np.sum(np.logical_and(self.wav < linewav + linewidth,self.wav > linewav - linewidth)) > 0)):
                        
                        # Add region
                        dellam = linewav*self.p['ContinuumRegion']/(2*C)
                        self.regions.append([np.max([linewav - dellam,self.wav[0]]), np.min([linewav + dellam,self.wav[-1]])])

                    # Else remove the line
                    else: 
                        for g in self.p['EmissionGroups']:
                            for s in g['Species']:
                                for l in s['Lines']:
                                    # Check if the same then remove
                                    if ((g['Name'] == group['Name']) and (s['Name'] == species['Name']) and (l['Wavelength'] == line['Wavelength'])):
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
        self.wav = self.wav[self.emission_region]
        self.flux = self.flux[self.emission_region]
        self.weight = self.weight[self.emission_region]
        self.sqrtweight = self.sqrtweight[self.emission_region]
        self.sigma = self.sigma[self.emission_region]

    # Boostrap the Flux
    def Boostrap(self):

        return np.random.normal(self.flux,self.sigma)