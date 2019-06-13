""" Find files and calculate cost functions.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import subprocess
import argparse
from smamp.tools import cd
from smamp.tools import check_existence

def calc_cost_function(path, lnrho=-9, sigma=0.8, output_file="cost.h5"):
	""" Call the horton script."""
	command = "horton-esp-cost.py"
	esp_file = " ../3_united_atom_structure/esp_ua.cube"
	output_file = ' ' + output_file
	density_file = " --wdens ../3_united_atom_structure/rho_ua.cube"
	density_file += ":{}:{}".format(lnrho, sigma)
	boundary_option = " --pbc 000 "
	sign_option = "--sign "  # changes sign
	overwrite = "--overwrite "

	# Build the command
	command += esp_file + output_file + density_file 
	command += boundary_option + sign_option + overwrite

	# Execute the command
	# print("Executing: \n{}".format(command))
	# print("In: {}".format(os.getcwd()))

	# output = subprocess.check_output(command, shell=True)
	with open('horton_stout.log', 'w') as logfile:
		kwargs = {"shell": True, "stdout": logfile, "stderr": subprocess.STDOUT}
		subprocess.Popen(command, **kwargs).communicate()
	# print('Calculation started, horton returned {}'.format(output.decode('ascii')))

def cmd_parser():
	parser = argparse.ArgumentParser(prog='',
					 description='Discover input files and calc cost functions.')
	
	parser.add_argument('--lnrho', metavar=-9,
        help='Value of the lnrho reference parameter for HORTON cost function construction.',
	default=-9)

	args = parser.parse_args()

	return args.lnrho


def main():
	""" Execute everything."""
	print('This is {}.'.format(__file__))
	print('Current working dir: {}'.format(os.getcwd()))

	lnrho =  cmd_parser()

	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk('.')):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folder to calculate in
		if 'horton_cost_function' in subdir:
			print('Moving to {}'.format(subdir))

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
				with cd(subdir):
					calc_cost_function(subdir, lnrho=lnrho)
	print('Done.')

if __name__ == '__main__':
	main()
