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
		if 'template' in subdir or 'exclude' in subdir:
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
	B_vector = np.array(f['cost']['B'])


	return A_matrix, B_vector

def collect_matrices(paths):
	"""
	Extract A and B from all cost functions.

	Arguments:
	paths: list of strings, pointing to the cost function files
	
	Returns:
	A_list: list of A-matrices
	B_list: list of B-vectors
	"""

	A_list = []
	B_list = []
	for cost_function_path in paths:
		A_matrix, B_vector = read_h5(cost_function_path)
		A_list += [A_matrix]
		B_list += [B_vector]

	assert len(A_list) == len(B_list)
	return A_list, B_list
	

def average(A_matrices, B_vectors):
	""" 
	Average over cost function matrices.

	The cost functions contain the three objects of the cost function: A, B, C
	A is a quadratic matrix (97x97), B a vector (d=97), and C is a constant.
	In the end, we are interested in the best-fit charges Q which are the solution to
	Q = A^-1 B

	Arguments:
	A_matrices: a list of NxN matrices
	B_vectors: a list of vectors with len N

	Returns:
	A: the average of all A_matrices.
	B: the average of all B_matrices.
	"""

	# Initialize empty
	A = A_matrices[0] * 0
	B = B_vectors[0] * 0

	# Average by adding all objects and dividing by their number
	for index in range(len(A_matrices)):
		A += A_matrices[index]
		B += B_vectors[index]

	# Divide
	number_snapshots = len(A_matrices)
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

	print('Data has been written to {}\n'.format(template_path))

def main():
	""" Run the script."""
	print('This is {}.'.format(__file__))
	# Create a folder for the averaged cost function
	chargepath = './horton_charges'
	if os.path.isdir(chargepath):
		shutil.rmtree(chargepath)
	os.mkdir(chargepath)

	# Find the locations of the cost function files
	cost_function_paths = find_cost()

	# Extract cost function As and Bs
	A_matrices, B_vectors = collect_matrices(cost_function_paths)

	# Average over all matrices & vectors
	average_A, average_B = average(A_matrices, B_vectors)

	# keep one HDF5 file as a template for writing into later
	shutil.copyfile(cost_function_paths[0], './horton_charges/costfunction_average.h5')

	# Export matrices to hdf5 
	export(average_A, average_B, template_path='./horton_charges/costfunction_average.h5')

	print('Done.')


if __name__ == '__main__':
	main()
