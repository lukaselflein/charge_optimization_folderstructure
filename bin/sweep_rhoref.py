"""Vary the rhoref parameter to find a sane value.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

from __future__ import print_function

import os
import shutil

from loop_cost_functions import calc_cost_function 
from smamp.tools import cd
from smamp.tools import check_existence

def iterate_lnrho(path_to_subdir):
	"""Vary the lnrho weighting parameter, create folder and execute."""
	sweep_dir = 'lnrho_sweep'
	if os.path.exists(sweep_dir):
		# print('Removing old dir')
		# shutil.rmtree(sweep_dir)
		pass
	else:
		print('making dir')
		os.mkdir(sweep_dir)
		print('dir made.')

	print('Iterating overt lnrho:')#, end=' ')
	for lnrho in range(-10, 0):
		output_name = os.path.join(sweep_dir, 'cost_{}.h5'.format(lnrho))
		print(output_name)#, end=' ')
		if os.path.exists(output_name):
			print('{} exists. Skipping ahead.')
			continue

		calc_cost_function(path_to_subdir, 
				   lnrho=lnrho, 
				   sigma=0.8, 
				   output_file=output_name)
	print()


def main():
	""" Execute everything."""
	print('This is {}.'.format(__file__))
	print('Current working dir: {}'.format(os.getcwd()))

	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk('.')):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folder to calculate in
		if 'horton_cost_function' in subdir:

			# Check if all neccessary files exist
			neccessary_files = ['../3_united_atom_structure/esp_ua.cube']
			neccessary_files += ['../3_united_atom_structure/rho_ua.cube']
			warning = check_existence(path=subdir, neccessary_files=neccessary_files)
	
			# If they don't exist, log this and move on
			if warning:
				with open('submissions.log', 'a') as logfile:
					logfile.write(warning + '\n')
				continue

			# Start the calculation only if all files exist
			elif warning is None:
				# print('Calculating: {}'.format(subdir))
				print('Moving to {}'.format(subdir))
				with cd(subdir):
					iterate_lnrho(subdir)	
	print('Done.')

if __name__ == '__main__':
	main()
