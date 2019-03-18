""" Search folderstructure for DFT output files, calculate bader charges.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import sys
import shutil
import argparse
import smamp
import subprocess

from smamp.tools import cd
from smamp.tools import find

def main():
	""" Execute everything."""
	print('This is {}.'.format(__file__))

	rho_paths = find(path='./', folder_keyword='dft', file_keyword='rho')

	for path in rho_paths:
		folder_path, file_name = os.path.split(path)
		command = 'bader -p atom_index '
		command += os.path.join(' ../2_dft_calculations/', file_name)
		topdir = os.path.split(folder_path)[0]
		bader_dir = os.path.join(topdir, '5_bader_charges')
		with cd(bader_dir):
			with open('bader.log', 'w') as logfile:
				kwargs = {"shell": True, "stdout": logfile, "stderr": subprocess.STDOUT}
				
				print('Executing bader in {}'.format(os.getcwd()))
				# Execute the command
				p = subprocess.Popen(command, **kwargs)
				# Wait for subprocess to finish
				_ = p.communicate()

				snapshot = '../0_initial_structure/snapshot.pdb'
				top = '../0_initial_structure/example.top'
				print('Extracting bader charges ...')
				smamp.extract_bader_charges.extract(snapshot, top)
				print('Extraction done.')

	print('Done.')

if __name__ == '__main__':
	main()
