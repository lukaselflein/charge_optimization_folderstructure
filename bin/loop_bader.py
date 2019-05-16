""" Search folderstructure for DFT output files, calculate bader charges.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import sys
import smamp
import subprocess

from smamp.tools import cd
from smamp.tools import find

def main():
	""" Execute everything."""
	print('This is {}.'.format(__file__))

	rho_paths = find(path='./', folder_keyword='dft', file_keyword='rho')

	for path in rho_paths:
		# Get the exact density file name
		folder_path, file_name = os.path.split(path)
		
		# Go to the bader directory to place the bader analysis output there	
		topdir = os.path.split(folder_path)[0]
		bader_dir = os.path.join(topdir, '5_bader_charges')
		with cd(bader_dir):
			print('Moving to {}'.format(bader_dir))

			# Find structure and topology files of the same snapshot
			snapshot_path = find(path='..', folder_keyword='initial', 
					     file_keyword='.pdb')[0]
			top_path = find(path='..', folder_keyword='initial', 
					file_keyword='.top')[0]

			# Write output to logfile	
			with open('bader.log', 'w') as logfile:
				# Assemble the shell command bader
				command = 'bader -p atom_index '
				command += os.path.join(' ../2_dft_calculations/', 
							file_name)
				kwargs = {"shell": True, "stdout": logfile, 
					  "stderr": subprocess.STDOUT}
				# Execute the shell command
				print('Running bader ...')
				#p = subprocess.Popen(command, **kwargs)

				# Wait for the shell command to finish
				#p.communicate()

				# Extract charges from the bader anaylsis output to .csv
				print('Bader done. Extracting bader charges ...')
				hyd_path = '../../../fitting_constraint_files/hydrogen_per_atom.csv'
				smamp.extract_bader_charges.extract(snapshot_path,
								    top_path,
								    hydrogen_path=hyd_path)
				print('Extraction done.')

	print('Done.')

if __name__ == '__main__':
	main()
