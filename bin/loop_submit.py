""" Search for inital united atom structures, convert them to all atom format.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import subprocess
from smamp.tools import cd

def check_existence(subdir):
	"""
	Check if all files needed for submission exist.
	"""
	with cd(subdir):
		# Submission files should be here
		print([f for s, d, f in os.walk('.')])
		neccessary_files = ['../1_all_atom_structure/ase_pdbH.traj']
		neccessary_files += ['gpaw_optimize_and_esp.py', 'submit_gpaw.sh']
		for f in neccessary_files:
			if not (os.path.islink(f) or os.path.isfile(f)):
				warning = 'File not found: {}'.format(f)
				print(warning)
				return warning
	return None

def submit(subdir):
	"""
	Submit the job description to the cluster queue.
	"""
	kwargs = {'stderr': subprocess.STDOUT, 'shell': True}
	output = subprocess.check_output('msub submit_gpaw.sh', **kwargs)
	print('Job submitted, queue response: {}'.format(output.decode('ascii')))

	with open('queue_response.log', 'ba+') as logfile:
		logfile.write(output)

def main():
	""" Execute everything."""
	print(__doc__)

	topdir = '.'

	# Crawl the directory structure
	for subdir, dirs, files in os.walk(topdir):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folders with submission scripts
		if 'dft_calculation' in subdir:
			print('Moving to {}'.format(subdir))

			# Check if all neccessary files exist
			warning = check_existence(subdir)
	
			# If they don't exist, log this and move on
			if warning:
				with open('submissions.log', 'a') as logfile:
					logfile.write(warning + '\n')
				# raise Warning(warning)
				continue

			# Only if all files exist, submit to the queue
			elif warning is None:
				print('Submitting {}'.format(subdir))
				with cd(subdir):
					submit(subdir)

if __name__ == '__main__':
	main()
