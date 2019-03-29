""" Search for 'all atoms' structures used in DFT, convert them to 'united atoms' format.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import smamp
from smamp.tools import cd
from smamp.tools import find


def convert(subdir):
	"""Convert united atom files to all atom format."""
	with cd(subdir):
		pdb_path = find(path='..', folder_keyword='initial', 
			        file_keyword='.pdb')[0]
		top_path = find(path='..', folder_keyword='initial', 
				file_keyword='.top')[0]

		for dft_type in ('esp', 'rho'):
			dft_path = find(path='..', folder_keyword='dft_calculations', 
					file_keyword='{}.cube'.format(dft_type))[0]

			out_file = dft_type + '_ua' + '.cube'
			kwargs = {'infile_pdb': pdb_path,
				  'infile_top': top_path,
				  'infile_cube': dft_path,
				  'outfile_cube': out_file}

			smamp.aa2ua_cube.aa2ua_cube(**kwargs)


def main():
	"""Execute everything."""
	print('This is {}.'.format(__file__))
	
	print('Current working dir: {}'.format(os.getcwd()))
	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk('.')):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folders with AA structures in them
		if 'united_atom_structure' in subdir:
			print('Moving to {}'.format(subdir))
			convert(subdir)


if __name__ == '__main__':
	main()
