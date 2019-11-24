""" Main Function """

# Packages
import sys
import copy
import argparse
import numpy as np
import multiprocessing as mp

# HQUAILS supporting files
import HQUAILS
import ConstructParams as CP

# Main Function
if __name__ == "__main__":

	## Parse Arguements to find Parameter File ##
	parser = argparse.ArgumentParser()
	parser.add_argument('Parameters', type=str, help='Path to parameters file')
	parser.add_argument('ObjectList', type=str, help='Path to object list with paths to spectra and their redshifts.')
	args = parser.parse_args()
	p = CP.construct(args.Parameters)
	## Parse Arguements to find Parameter File ##

	## Verify Emission Line Dictionary ##
	if not CP.verify(p['EmissionLines']):
		print('Unable to verify emission line dictionary, exiting.')
		sys.exit(1)
	## Verify Emission Line Dictionary ##

	## Assemble Objects
	objects = np.genfromtxt(args.ObjectList,delimiter=',',dtype='U100,f8',names=['File','z'])
	## Assemble Objects

	## Run HQUAILS ##
	if p['NPool'] == 1: # Single Thread
		for o in objects: HQUAILS.HQUAILS(copy.deepcopy(p),o['File'],o['z'])
	else: # Multithread
		pools = mp.Pool(processes=p['NPool'])
		inputs = [(copy.deepcopy(p),o['File'],o['z']) for o in objects]
		pools.starmap(HQUAILS.HQUAILS, inputs)
	## Run HQUAILS ##
