#!/usr/bin/env python

import os
from astropy.io import fits
from CustomModels import SSPContinuum
from SpectrumClass import Spectrum
from ConstructParams import construct
from FittingModel import FitModel
import numpy as np

import matplotlib.pyplot as plt

params = construct('Example/ExPARAMS.json')
wav = 10**fits.getdata('Example/Spectra/spec-0295-51585-0596.fits')['loglam']


spectrum = Spectrum('Example/Spectra/spec-0295-51585-0596.fits',0.138345125,params)
print(spectrum.emission_region.sum())
# continuum = SSPContinuum(spectrum)
# region = np.invert(region)
# continuum.set_region(region)

# fit_model = FitModel(spectrum,continuum,region)
# fit_model.set_region(np.ones(spectrum.wav.shape,dtype=bool))


# # plt.plot(spectrum.wav,spectrum.flux,ds='steps-mid',label='Data',color='k',lw=1)
# # plt.plot(wav,fit_model(wav),ds='steps-mid',label='Model',color='r',lw=1)

# plt.plot(spectrum.wav,spectrum.flux-fit_model(wav),ds='steps-mid',label='Residuals',color='r',lw=0.5)
# plt.plot(spectrum.wav,spectrum.sigma**2,ds='steps-mid',label='Variance',color='k',lw=0.5)
# # plt.xlim([7000,8000])
# # plt.ylim([0,None])
# plt.legend()
# plt.tight_layout()
# plt.savefig('test.pdf')