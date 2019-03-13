""" Search for 'all atoms' structures used in DFT, convert them to 'united atoms' format.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import sys
import shutil
from smamp.tools import cd
import smamp
import warnings

def convert(subdir):
	"""
	Converts united atom files to all atom format.
	"""
	with cd(subdir):
		# The original UA file should be in ../0_initial_structure
		init_path = '../0_initial_structure'
		# The DFT files should be in ../2_dft_calculations
		dft_path = '../2_dft_calculations'
		for name in ('.pdb', '.top'):
			for subdir, dirs, files in os.walk(init_path):
				# Check if at least one file exists
				if sum([name in f for f in files]) < 1:
					print(os.getcwd())
					raise RuntimeError('No {} file found.'.format(name))
				# No more than one file must exists for uniqueness
				if sum([name in f for f in files]) > 2:
					print(os.getcwd())
					raise RuntimeError('Multiple {} files found.'.format(name))

		name = '.cube'
		for subdir, dirs, files in os.walk(dft_path):
			# Check if at least one file exists
			if sum([name in f for f in files]) < 2:
				print(os.getcwd())
				print(dirs, files)
				raise RuntimeError('Not enough {} files found.'.format(name))
			# No more than one file must exists for uniqueness
			if sum([name in f for f in files]) > 3:
				print(os.getcwd())
				print(dirs, files)
				raise RuntimeError('Too many {} files found.'.format(name))

		for dft_file in ('esp.cube', 'rho.cube'):
			out_file = dft_file[:3] + '_ua' + dft_file[3:]
			kwargs = {'infile_pdb': os.path.join(init_path, 'snapshot.pdb'),
				  'infile_top': os.path.join(init_path, 'example.top'),
				  'infile_cube': os.path.join(dft_path, dft_file),
				  'outfile_cube': out_file}

			smamp.aa2ua_cube.aa2ua_cube(**kwargs)

def main():
	"""Execute everything."""
	print('This is {}.'.format(__file__)
	# Save the working dir
	topdir = os.getcwd()
	print('Current working dir: {}'.format(topdir))
	
	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk(topdir)):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folders with AA structures in them
		if 'united_atom_structure' in subdir:
			print('Moving to {}'.format(subdir))
			convert(subdir)


if __name__ == '__main__':
	main()
