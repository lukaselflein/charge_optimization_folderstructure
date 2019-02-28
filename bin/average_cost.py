""" Average Cost Functions for Horton to determine Charges for Molecular Dynamics.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import numpy as np
import h5py
import shutil
import warnings
import os
from smamp.tools import cd


def find_cost(path='.', cost_function_filename='cost.h5'):
	"""
	Find all cost functions in folder structure.

	We have multiple snapshots of a molecule, with corresponding cost functions.
	This function explores the folderstructure, and returns all paths leading to cost functions.

	Arguments:
	path: path to the folder to search in

	Returns:
	list of strings, representing paths to cost functions
	"""

	cost_function_paths = []

	# Crawl the directory structure
	for subdir, dirs, files in os.walk(path):

		# Exclude template folders from search
		if 'template' in subdir:
			continue

		# Select the folder with cost functions:
		if 'horton_cost_function' in subdir:
			# The cost file should be in:
			cost_path = os.path.join(subdir, cost_function_filename)
			if os.path.isfile(cost_path):
				# Add the cost function to our collection
				cost_function_paths += [cost_path]

			# Maybe the cost function is missing. Print to screen
			else:
				print('\nWarning: No cost file found in: "{}"'.format(subdir))
				print('Filename was assumed to be "{}"'.format(cost_function_filename))
	return cost_function_paths

def read_h5(path, verbose=True):
	"""
	Import cost functions, convert to numpy arrays.

	Argumentis:
	path: the path to the cost function file

	Returns:
	A_matrices: a dictionary of matrices, indexed by their timestep.
	B_vectors: a dictionary of vectors, indexed by their timestep.
	"""


	# Extract the values for each timestep
	if verbose:
		print('Loading data: {}'.format(path))

	# load the objects (read-only) from HDF5 file
	f = h5py.File(path, 'r')
	# Extract the A matrix
	A_matrix = np.array(f['cost']['A'])
	# Extract the B vector
	B_matrix = np.array(f['cost']['B'])


	return A_matrix, B_vector

def collect_matrices():
	pass

def average(A_matrices, B_vectors, timesteps):
	""" 
	Average over cost function matrices.

	The cost functions contain the three objects of the cost function: A, B, C
	A is a quadratic matrix (97x97), B a vector (d=97), and C is a constant.
	In the end, we are interested in the best-fit charges Q which are the solution to
	Q = A^-1 B

	Arguments:
	A_matrices: a dictionary of NxN matrices, indexed with their timestep.
	B_vectors: a dictionary of vectors with len N, indexed with their timestep.
	timesteps: a list or tuple of timesteps.

	Returns:
	A: the average of all A_matrices.
	B: the average of all B_matrices.
	"""

	# Initialize empty
	time = list(B_vectors.keys())[0]
	A = A_matrices[time] * 0
	B = B_vectors[time] * 0

	# Average by adding all objects and dividing by their number
	for timestep in timesteps:
		A += A_matrices[timestep]
		B += B_vectors[timestep]

	# Divide
	number_snapshots = len(timesteps)
	A /= number_snapshots
	B /= number_snapshots

	return A, B

def export(A, B, template_path='./average_cost.h5'):
	"""
	Export&save numpy-matrices to HDF5 objects.
	
	Arguments:
	A: Averaged A-matrix.
	B: Averaged B-vector.

	Keyword Arguments:
	template_path: the path to the template file A and B are written into.
	"""
	# Open the template file
	f = h5py.File(template_path, 'r+')
	# Load the template A matrix
	A_old = f['cost/A']
	# Assign the averaged A
	A_old[...] = A
	# Do the same for the B-vectors
	B_old = f['cost/B']
	B_old[...] = B
	# Save changes
	f.close() 

	# Make sure that the changes were written into the template
	f = h5py.File(template_path, 'r')
	assert np.allclose(f['cost/A'][()], A)
	assert np.allclose(f['cost/B'][()], B)

	print('\nData has been written to {}:'.format(template_path))

def main():
	"""
	Run the script.
	"""
	cost_function_paths = find_cost()	

	# keep one HDF5 file as a template for writing into later
	shutil.copyfile(cost_function_paths[0], './average_cost.h5')

	
	# WORK_DIR = '.'
	# TIMESTEPS = [str(time) for time in range(100, 1100, 100)] 
	# A_matrices, B_vectors = read_h5(work_dir=WORK_DIR, timesteps=TIMESTEPS)

	# print('Calculating averages')
	# A, B = average(A_matrices, B_vectors, timesteps=TIMESTEPS)
	
	# export(A, B)
	print('Done.')

if __name__ == '__main__':
	main()
