""" Search for inital united atom structures, convert them to all atom format.
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

def calc_charges(subdir, qtot=0):
	""" 
	Calculates HORTON charges for one cost function. 

	A call to `fitESPconstrained.py` is assembled via path-strings, and executed.
	"""
	# Specify the filepaths
	binary = 'python ../../../bin/fitESPconstrained.py '
	target_cost = 'cost.h5 '
	snapshot_path = '../0_initial_structure/snapshot.pdb '
	topology_path =  '../0_initial_structure/example.top '
	constraint_path = '../../../fitting_constraint_files/'
	constraint_files = constraint_path + 'atoms_in_charge_group.csv ' 
	constraint_files += constraint_path + 'charge_group_total_charge.csv '
	constraint_files += constraint_path + 'atoms_of_same_charge.csv '
	output_files = 'fitted_point_charges.csv '
	# outputfiles += 'test.top '
	# Specify options
	charge_option = '--qtot {} '.format(qtot)
	verbosity = '--verbose '

	# Assemble the command
	command = binary + target_cost + snapshot_path + topology_path
	command += constraint_files + output_files + charge_option + verbosity

	# Log the command output
	with open('charge_fitting_out.log', 'w') as logfile:
		kwargs = {"shell": True, "stdout": logfile, "stderr": subprocess.STDOUT}
		# Execute the command
		p = subprocess.Popen(command, **kwargs).communicate()
	print(p)

def main():
	""" Execute everything."""

	# TODO: make this a command line argument
	qtot=0

	# Save the working dir
	topdir = os.getcwd()
	print('Current working dir: {}'.format(topdir))
	
	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk(topdir)):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folders with cost function
		if 'horton_cost_function' in subdir:
			print('Moving to {}'.format(subdir))
			with cd(subdir):
				calc_charges(subdir, qtot=qtot)


if __name__ == '__main__':
	main()
