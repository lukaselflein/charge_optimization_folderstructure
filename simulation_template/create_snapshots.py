""" Create snapshots from trajectory, initialize folder structure 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import shutil
import argparse
import os

def cmd_parser():
	"""
	Read Command line arguments
	"""
	parser = argparse.ArgumentParser(prog='',
					 description='Create snapshots from trajectory in folder structure')
	parser.add_argument('-wd',
        help='The top-level working directory containing the snapshot subfolders',
	default='./', metavar='./')

	parser.add_argument('-n', metavar='10',
        help='The number of snapshots to be created',
	default='10', type=int)

	parser.add_argument('-s', metavar='100',
        help='The timestamp of the first snapshot',
	default='100', type=int)

	parser.add_argument('-d', metavar='100',
        help='The difference in time between snapshots',
	default='100', type=int)

	args = parser.parse_args()

	return args.wd, args.n, args.s, args.d

def main():
	"""
	Run the script.
	"""

	# Read Command-Line Arguments
	wd, number_snapshots, start_time, delta_time = cmd_parser()

	# Calculate the corresponding timesteps
	end_time = start_time + delta_time * (number_snapshots)
	timesteps = range(start_time, end_time, delta_time)
	print('Creating directories for timesteps {}'.format(timesteps))

	# Set up directories for each timestep
	for time in timesteps:
		path = wd + str(time) + '_ps_snapshot'
		
		# Make sure no old snapshot dir exists - if it does, remove it.
		if os.path.isdir(path):
			print('{} alread exists. Deleting it ...'.format(path))
			shutil.rmtree(path)

		# Copy the template
		shutil.copytree(wd + 'snapshot_template', path, symlinks=True)
		print('{} initialized.'.format(path))

	print('Done.')

if __name__ == '__main__':
	main()
