""" Search for inital united atom structures, convert them to all atom format.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import sys
import shutil
import argparse
import smamp

from smamp.tools import cd
from smamp.tools import read_atom_numbers


def convert(subdir, hydrogen_per_atom):
	"""Converts united atom files to all atom format."""
	with cd(subdir):
		# The UA file should be in ../0_initial_structure
		path = '../0_initial_structure'
		for subdir, dirs, files in os.walk(path):
			# Check if at least one file exists
			if sum(['snapshot' in f for f in files]) < 1:
				print(os.getcwd())
				raise RuntimeError('No snapshot.pdb file found.')
			# No more than one file must exists for uniqueness
			if sum(['snapshot' in f for f in files]) > 2:
				print(os.getcwd())
				raise RuntimeError('Multiple snapshot.pdb files found.')

		# Convert the united atoms file
		smamp.convert_UA_to_AA.main(implicitHbondingPartners=hydrogen_per_atom)

def cmd_parser():
	"""Read Command line arguments."""
	parser = argparse.ArgumentParser(prog='',
			description=__doc__.splitlines()[0])

	args = parser.parse_args()
	

def main():
	""" Execute everything."""
	cmd_parser()
	print('This is {}.'.format(__file__))

	# Read number of hydrogen per atom from configuration file table
	hydrogen_per_atom = read_atom_numbers()

	print('Current working dir: {}'.format(os.getcwd()))
	
	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk('.')):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folders with AA structures in them
		if 'all_atom_structure' in subdir:
			print('Moving to {}'.format(subdir))
			convert(subdir, hydrogen_per_atom)


if __name__ == '__main__':
	main()
