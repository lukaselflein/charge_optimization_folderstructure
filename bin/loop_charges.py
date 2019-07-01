""" Calculate HORTON charges for all .cube files in folderstructure.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import sys
import shutil
import smamp
import warnings
import subprocess

from smamp.tools import cd
from smamp.tools import read_total_charge
from smamp.tools import find


def calc_charges(subdir, qtot=0, pdb_path='../0_initial_structure/snapshot.pdb', top_path='../0_initial_structure/example.top'):
	"""Calculates HORTON charges for one cost function. 

	A call to `fitESPconstrained.py` is assembled via path-strings, and executed.
	"""
	# Specify the filepaths
	binary = 'python ../../../bin/fitESPconstrained.py '
	target_cost = 'cost.h5 '

	snapshot_path = pdb_path + ' '
	topology_path =  top_path + ' '
	constraint_path = '../../../fitting_constraint_files/'
	constraint_files = constraint_path + 'atoms_in_charge_group.csv ' 
	constraint_files += constraint_path + 'charge_group_total_charge.csv '
	constraint_files += constraint_path + 'atoms_of_same_charge.csv '
	output_files = 'fitted_point_charges.csv '
	hydrogen_file = ' -i ' + constraint_path + 'hydrogen_per_atom.csv '
	# outputfiles += 'test.top '
	# Specify options
	charge_option = '--qtot {} '.format(qtot)
	verbosity = '--verbose '

	# Assemble the command
	command = binary + target_cost + snapshot_path + topology_path
	command += constraint_files + output_files + charge_option 
	command += hydrogen_file + verbosity

	# Log the command output
	with open('charge_fitting_out.log', 'w') as logfile:
		kwargs = {"shell": True, "stdout": logfile, "stderr": subprocess.STDOUT}
		# Execute the command
		p = subprocess.Popen(command, **kwargs).communicate()

def main():
	""" Execute everything."""
	print('This is {}.'.format(__file__))
	# TODO: make this a command line argument

	# Save the working dir
	topdir = os.getcwd()
	print('Current working dir: {}'.format(topdir))

	# Read the file containing the info on the total charge of the system
	qtot = read_total_charge(path='../fitting_constraint_files/total_charge.csv')
	
	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk('.')):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir or 'sweep' in subdir:
			continue

		# Select the folders with cost function
		if 'horton_cost_function' in subdir:
			print('Moving to {}'.format(subdir))
			with cd(subdir):
				# Search for the .pdb and .top file in the folders above
				input_path = '..'
				pdb_path = find(input_path, folder_keyword='initial', file_keyword='.pdb', nr_occ=1, exclude_kw=['template', 'exclude', 'sweep'])[0]
				top_path = find(input_path, folder_keyword='initial', file_keyword='.top', nr_occ=1, exclude_kw=['template', 'exclude', 'sweep'])[0]
	
				# Do the calculation
				calc_charges(subdir, qtot=qtot, pdb_path=pdb_path, top_path=top_path)


if __name__ == '__main__':
	main()
