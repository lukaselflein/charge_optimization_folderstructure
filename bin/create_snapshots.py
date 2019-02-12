""" Create snapshots from trajectory, initialize folder structure 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import shutil
import argparse
import os
import subprocess

def cmd_parser():
	"""
	Read Command line arguments
	"""
	parser = argparse.ArgumentParser(prog='',
					 description='Create snapshots from trajectory in folder structure')
	parser.add_argument('-wd',
        help='The top-level working directory containing the snapshot subfolders',
	default='./', metavar='./')

	parser.add_argument('-trj',
        help='The path to the trajectory file',
	default='./traject.xtc', metavar='./traject.xtc')

	parser.add_argument('-tpr',
        help='The path to the topolgy file',
	default='./topology.tpr', metavar='./topo.tpr')

	parser.add_argument('-e', metavar='10',
        help='The timestamp of the last snapshot', 
	default='50', type=int)

	parser.add_argument('-s', metavar='100',
        help='The timestamp of the first snapshot',
	default='0', type=int)

	parser.add_argument('-d', metavar='100',
        help='The difference in time between snapshots',
	default='10', type=int)

	args = parser.parse_args()

	return args.wd, args.e, args.s, args.d, args.trj, args.tpr

def create_directories(working_dir, timesteps):
	"""
	Set up directories for each timestep
	"""
	for time in timesteps:
		path = working_dir + str(time) + '_ps_snapshot'
		
		# Make sure no old snapshot dir of same name exists - if it does, remove it.
		if os.path.isdir(path):
			print('{} alread exists. Deleting it ...'.format(path))
			shutil.rmtree(path)

		# Copy the template
		shutil.copytree(working_dir + 'snapshot_template', path, symlinks=True)
		print('{} initialized.'.format(path))

	return path

def extract_snapshots(trajectory_file, tpr_file, path, start_time, end_time, delta_time):
	print(trajectory_file)
	snapshot_name = 'snapshot.pdb'
	command = 'gmx trjconv -f {xtc} -s {tpr}'.format(xtc=trajectory_file, tpr=tpr_file)
	command += '-o {pdb} -pbc mol -b {start} -e {end} -dt {dt} -sep'.format(pdb=snapshot_name, start=start_time, end=end_time, dt=delta_time)
	try:
		subprocess.call(command)
	except:
		raise RuntimeError('"{}" failed, no snapshots were extracted.'.format(command))

	first_name = 'snapshot0.pdb'
	if os.path.isfile(os.join(path, first_name)):
		print(first_name)
	
	print('Snapshots successfully extracted.')

def main():
	"""
	Run the script.
	"""

	# Read Command-Line Arguments
	working_dir, end_time, start_time, delta_time, trajectory_file, tpr_file = cmd_parser()

	# Calculate the corresponding timesteps
	timesteps = range(start_time, end_time, delta_time)
	print('Creating directories for timesteps {}'.format(timesteps))
	path = create_directories(working_dir, timesteps)
	extract_snapshots(trajectory_file=trajectory_file, tpr_file=tpr_file, 
			  path=path, start_time=start_time, delta_time=delta_time, 
			  end_time=end_time)

	print('Done.')

if __name__ == '__main__':
	main()
