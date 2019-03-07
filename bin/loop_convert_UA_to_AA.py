""" Search for inital united atom structures, convert them to all atom format.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import sys
import shutil

# Hacky relative import to avoid having to install a package
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'bin'))
import convert_UA_to_AA

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def convert(subdir):
	"""
	Converts united atom files to all atom format.
	"""
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
		convert_UA_to_AA.main()
	

def main():
	"""
	Execute everything.
	"""

	# Save the working dir
	topdir = os.getcwd()
	print('Current working dir: {}'.format(topdir))
	
	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk(topdir)):

		# Exclude template folders from search
		if 'template' in subdir:
			continue

		# Select the folders with AA structures in them
		if 'all_atom_structure' in subdir:
			print('Moving to {}'.format(subdir))
			convert(subdir)


if __name__ == '__main__':
	main()
