""" Average Cost Functions for Horton to determine Charges for Molecular Dynamics.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import numpy as np
import h5py
import shutil
import warnings
import os
from smamp.tools import cd
from smamp.tools import find
from average_cost import read_h5, collect_matrices, average, export


def main():
	""" Run the script."""
	print('This is {}.'.format(__file__))
	# Create a folder for the averaged cost function
	chargepath = './horton_charges/sweep_rhoref'
	if os.path.isdir(chargepath):
		pass
		# shutil.rmtree(chargepath)
	else:
		os.makedirs(chargepath)

	cost_paths = find(path='.', folder_keyword='4_horton_cost_function/lnrho_sweep', 
			  file_keyword='cost',
			  nr_occ=None)
	lnrho_range = []
	sigma_range = []
	for charge_file in cost_paths:
		# print(charge_file)
		# Parse parameters from filename
		lnrho, sigma = charge_file[-15:-3].split('_')[-2:]
		lnrho_range += [lnrho]
		sigma_range += [sigma]

	lnrho_range = set(lnrho_range)
	sigma_range = set(sigma_range)

	for lnrho in lnrho_range:
		for sigma in sigma_range:
			# print('lnrho = {}, sigma = {}'.format(lnrho, sigma))
			filename = 'cost_{}_{}.h5'.format(lnrho, sigma)
			cost_function_paths = find(path='.', folder_keyword='lnrho_sweep', 
						   file_keyword=filename, nr_occ=None)

			# Extract cost function As and Bs
			A_matrices, B_vectors = collect_matrices(cost_function_paths)

			# Average over all matrices & vectors
			average_A, average_B = average(A_matrices, B_vectors)

			# keep one HDF5 file as a template for writing into later
			out_name =  'costfunction_average_{}_{}.h5'.format(lnrho, sigma)
			template_path = os.path.join(chargepath, out_name) 
			shutil.copyfile(cost_function_paths[0], template_path)

			# Export matrices to hdf5 
			export(average_A, average_B, template_path=template_path)

	print('Done.')


if __name__ == '__main__':
	main()
