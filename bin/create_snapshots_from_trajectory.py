""" Create snapshots from trajectory, initialize folder structure 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import shutil
import argparse
import os
import subprocess

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
		shutil.copytree(working_dir + '.template_snapshot', path, symlinks=True)
		print('{} initialized.'.format(path))

	return path

def extract_snapshots(trajectory_file, tpr_file, top_file, working_dir, 
                      start_time, end_time, delta_time):
	"""
	Extract snapshots via GROMACS and move them to their snapshot-subfolders.

	We use the trajectory of the molecule in solution to extract structures at different time steps.
	The exact time steps are given by start, end and delta time.

	Args:
	trajectory_file: Path to the GROMACS .xtc trajectory file
	tpr_file: Path to the topology file
	working_dir: Path to the snapshots-subfolder
	"""

	# Compose the command to GROMACS
	snapshot_name = 'snapshot.pdb'
	# Select the "other" group containing the molecule, excluding the solution:
	command = 'echo 1 | '
	# Specify input files
	command += 'gmx trjconv -f {xtc} -s {tpr}'.format(xtc=trajectory_file, tpr=tpr_file)
	# Specify output and timesteps
	command += ' -o {pdb} -pbc mol -b {start} -e {end} -dt {dt} -sep'.format(pdb=snapshot_name, 
										start=start_time,
										end=end_time,
										dt=delta_time)

	print('\nRunning GROMACS trjconv, please wait ...')

	verbose = True
	if verbose == True:
		try:
			subprocess.check_call([command], shell=True)
		except:
			raise RuntimeError('I got an error from: \n "{}"\
					   \n Did you load GROMACS?'.format(command))

	else:
		# This call catches the stderr output without printing it
		pipe = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, 
					stderr=subprocess.PIPE)
		result = pipe.communicate()
		with open('gmx_trjconv.log', 'wb') as logfile:
			for line in result:
				logfile.write(line)
		print('GROMACS finished successfully.')

	# Check if output exists, then move to respective folders
	index = 0  # This is for GROMACS naming convention in [0..n]
	for time in range(start_time, end_time + 1, delta_time):
	
		# Define names
		snapshot_name = 'snapshot' + str(index) + '.pdb' 
		snapshot_path = os.path.join(working_dir, snapshot_name)
		target_folder = os.path.join(working_dir, str(time) + '_ps_snapshot', 
					     '0_initial_structure')

		# Check if snapshots was actually created
		if not os.path.isfile(snapshot_path):
			raise RuntimeError('Snapshot {} is missing'.format(snapshot_path))
		
		# Check if folder exists, e.g. 600_ps_snapshot
		if not os.path.isdir(target_folder):
			raise RuntimeError('Folder {} is missing'.format(target_folder))

		# Move snapshot to its own folder
		shutil.move(snapshot_path, os.path.join(target_folder, 'snapshot.pdb'))

		# Copy the .top to the subfolders
		top_name = os.path.split(top_file)[-1]
		shutil.copyfile(top_file, os.path.join(target_folder, top_name))
		
		index += 1  # for Gromacs naming
	
	print('Snapshots successfully extracted.')

def cmd_parser():
	"""
	Read Command line arguments
	"""
	parser = argparse.ArgumentParser(prog='',
					 description='Create snapshots from trajectory in folder structure')
	parser.add_argument('-wd',
        help='The top-level working directory containing the snapshot subfolders',
	default='./', metavar='./')

	parser.add_argument('-xtc',
        help='The path to the trajectory file',
	default='./simulation/example.xtc', metavar='./traject.xtc')

	parser.add_argument('-tpr',
        help='The path to the topolgy file',
	default='./simulation/example.tpr', metavar='./topo.tpr')

	parser.add_argument('-top',
        help='The path to the .top topolgy file',
	default='./simulation/example.top', metavar='./topo.top')

	parser.add_argument('-e', metavar='1000',
        help='The timestamp of the last snapshot', 
	default='1000', type=int)

	parser.add_argument('-s', metavar='100',
        help='The timestamp of the first snapshot',
	default='200', type=int)

	parser.add_argument('-d', metavar='100',
        help='The difference in time between snapshots',
	default='100', type=int)

	args = parser.parse_args()

	return args.wd, args.e, args.s, args.d, args.xtc, args.tpr, args.top

def main():
	"""Run the script."""
	print('This is {}.'.format(__file__)

	# Read Command-Line Arguments
	working_dir, end_time, start_time, delta_time, trajectory_file, tpr_file, top_file = cmd_parser()

	# Calculate the corresponding timesteps
	timesteps = range(start_time, end_time + delta_time, delta_time)
	print('Creating directories for timesteps {}'.format(timesteps))
	create_directories(working_dir, timesteps)
	extract_snapshots(trajectory_file=trajectory_file, tpr_file=tpr_file, top_file=top_file,
			  working_dir=working_dir, start_time=start_time, delta_time=delta_time, 
			  end_time=end_time)


	print('\nAll done.')

if __name__ == '__main__':
	main()
