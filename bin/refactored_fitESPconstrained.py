""" Fit point charges to a HORTON costfunction under constraints.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
Based on legacy code by Johannes Hormann
"""

import argparse
import h5py
import warnings
import ase.io

import parmed as pmd
import numpy as np
import pandas as pd

from smamp.insertHbyList import insertHbyList
from smamp.tools import read_atom_numbers

def create_structure(infile_pdb, infile_top, strip_string=':SOL,CL', implicitHbondingPartners = None):
        
    if implicitHbondingPartners is None:
        implicitHbondingPartners = read_atom_numbers()
        
    ua_ase_struct = ase.io.read(infile_pdb)
    ua_pmd_struct = pmd.load_file(infile_pdb)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ua_pmd_top = pmd.gromacs.GromacsTopologyFile(infile_top,parametrize=False)

    ua_pmd_top.strip(strip_string)
        #strip water and electrolyte from system (if not yet done in .top)
    ua_pmd_top.box = ua_pmd_struct.box # Needed because .pdb contains box info
    ua_pmd_top.positions = ua_pmd_struct.positions

    ua_names = [ a.name for a in ua_pmd_top.atoms ]
    ua_residues = [ a.residue.name for a in ua_pmd_top.atoms ]

    aa_ase_struct, aa_pmd_struct, aa_names, aa_residues = \
        insertHbyList(ua_ase_struct,ua_pmd_top,
        implicitHbondingPartners,1.0)

    ua_count = len(ua_ase_struct)     # united atoms structure
    aa_count = len(aa_ase_struct) # all atoms structure

    ua_ase_index = np.arange(ua_count)
    aa_ase_index = np.arange(aa_count)

    aa_atom_residue_list = list(zip(aa_names,aa_residues))
    aa_ase_index = range(aa_count)
    aa_ase2pmd = dict(zip(aa_ase_index,aa_atom_residue_list))
    aa_pmd2ase = dict(zip(aa_atom_residue_list,aa_ase_index))

    ua_atom_residue_list = list(zip(ua_names,ua_residues))
    ua_ase_index = range(ua_count)
    ua_ase2pmd = dict(zip(ua_ase_index,ua_atom_residue_list))
    ua_pmd2ase = dict(zip(ua_atom_residue_list,ua_ase_index))

    # TODO: distinction for ua and aa fitting:
    pmd_struct = ua_pmd_struct
    pmd_top = ua_pmd_top
    ase2pmd = ua_ase2pmd
    pmd2ase = ua_pmd2ase
    return pmd_struct, pmd_top, ase2pmd

def constrained_minimize(A, B, D, Q):

	# Default to zero total charge constraint
	if D is None and Q is None:
		Q = np.array([0])
		D = np.ones(B.shape[0])
		

	# Minimize the cost
	# cost = xAx - 2Bx - C 

	# Cast everything to arrays
	A = np.atleast_2d(A)
	D = np.atleast_2d(D)
	B = np.atleast_1d(B)
	Q = np.atleast_1d(Q) 

	# For old versions of numpy, block is not available. Fallback to bmat:
	if  float(np.version.version[2:]) < 13:
		stack = np.bmat
	else:
		stack = np.block

	# Stack the HORTON matrices with the constraints
	# Unconstrained:
	# A x - B = 0
	#
	# With constraints:
	# A D  x  =  B  
	# D 0  l     Q
	zeros = np.zeros((Q.shape[0], Q.shape[0]))
	A_con = stack([[A, D.T],
		       [D, zeros]])
	B_con = stack([B, Q]).T

	x = np.linalg.solve(A_con, B_con)

	charges = x[:len(B)]
	langrange_forces = x[len(B):]
	return charges, langrange_forces


def parse_charge_groups(file_name, ase2pmd):
	# first we read in the textfile
	df = pd.read_csv(file_name, sep=',', header=None,
			 comment='#', names=['atom','cg'])

	# Atoms appear in multiple charge groups.
	# In the end, we want something like
	# {cg1: [1, 5, 8]}
	charge_groups = {}
	for atom in df.atom:
		# cg is the charge group of the current atom
		cg = df.loc[df.atom == atom].cg.values[0]
		if not cg in charge_groups.keys():
			charge_groups[cg] = []        

		# ase2pmd is formatted like
		# 0: ('CE1', 'terB')
		for ase_index, atom_residuum in ase2pmd.items():
			# If the atom names match, pick the ase index            
			if atom in atom_residuum:
				charge_groups[cg] += [ase_index]
	# Sort everything                
	for ase_index in charge_groups.keys():
		charge_groups[ase_index].sort()

	return charge_groups

def parse_group_charges(file_name):
	group_q = pd.read_csv(file_name, sep=',', header=None, comment='#',
	                      names=['charge'], index_col=0)
	group_q.charge = group_q.charge.astype(int)    
	return group_q

def parse_symmetry(path):
	print(path)


def make_group_constraints(charge_groups, group_q, n_atoms):
	# Initialize empty arrays
	D_matrix = np.zeros((len(charge_groups), n_atoms), dtype=int)
	Q_vector = np.zeros(len(charge_groups), dtype=int)

	# Fill in constraint values for every charge group
	for ase_index in charge_groups.keys():
		cg = charge_groups[ase_index]
		# Note: Charge groups are [1, 2, ...], np indices are [0, 1, ..]
		# 1 means that the sum of q_i in the charge group is unweighted
		D_matrix[ase_index - 1, cg] = 1

		# Now we need to specify the total charge of the group in a vector
		total_group_charge = group_q.loc[ase_index]
		Q_vector[ase_index - 1] = total_group_charge
	return D_matrix, Q_vector


def stack_constraints(group_matrix, group_q, symmetry_matrix, symmetry_q):

	if not None in (group_matrix, group_q, symmetry_matrix, symmetry_q):
		constraint_matrix = np.vstack(group_matrix, symmetry_matrix)
		constraint_q = np.vstack(group_q, symmetry_q)
		return constraint_matrix, constraint_q

	if group_matrix is None and symmetry_matrix is not None:
		return symmetry_matrix, symmetry_q

	if group_matrix is not None and symmetry_matrix is None:
		return group_matrix, group_q

	else:
		return None, None


def get_constraints(args, ase2pmd):
	'''Read provided constraint files and convert them into matrix form.'''
	charge_group_file = args.charge_groups
	charge_group_charges_file = args.charge_group_charges
	symmetry_file = args.symmetry_file
	
	if charge_group_file is not None:
		if charge_group_charges_file is None:
			err = 'Charge groups defined: {}'.format(charge_group_file)
			err += '\n But no total charges were defined.'
			raise ValueError(err)

		charge_groups = parse_charge_groups(charge_group_file, ase2pmd)
		group_q = parse_group_charges(charge_group_charges_file)
		group_matrix, group_q = make_group_constraints(charge_groups, group_q, n_atoms)

	else:
		group_matrix, group_q = None, None

	if symmetry_file is not None:
		symmetry = parse_symmetry(symmetry_file)
		symmetry_matrix, symmetry_q = make_symmetry_constraints(symmetry_file)
	else:
		symmetry_matrix, symmetry_q = None, None

	# Remove redunant constraints
	constraint_matrix, constraint_q = stack_constraints(group_matrix, group_q,
				                            symmetry_matrix, symmetry_q)
	
	return constraint_matrix, constraint_q


def read_horton_cost_function(file_name):
    cost_function = h5py.File(file_name)
    A = cost_function['cost']['A'][()]
    B = cost_function['cost']['B'][()]
    return A, B


def write_charges(charges):
	pass


def write_forces(f):
	pass


def get_ase_indices():
	pass

def parse_command_line():
	"""Read file locations from command line interface."""
	parser = argparse.ArgumentParser(prog='esp-fit-constrained.py',
				         description='Estimate charges from a HORTON ESP'
						     'cost function under constraints.')
	
	parser.add_argument('-hor', '--horton_cost_function',
		help='The location of the HORTON cost function file.',
		required=True, metavar='cost.h5')

	parser.add_argument('-p', '--pdb_infile',
		help='The location of the atomic structure file',
		required=True, metavar='snapshot.pdb')

	parser.add_argument('-t', '--top_infile',
		help='The location of the topolgy file',
		required=True, metavar='topol.top')

	parser.add_argument('-g', '--charge_groups',
		help='The location of the charge group constraints .csv file.',
		metavar='atoms_in_charge_group.csv', default=None)

	parser.add_argument('-c', '--charge_group_charges',
		help='The location of the charge group total charges .csv file.',
		metavar='charge_group_total_charge.csv', default=None)

	parser.add_argument('-s', '--symmetry_file',
		help='The location of the symmetry constraints file.',
		metavar='atoms_of_same_charge.csv', default=None)

	parser.add_argument('-o', '--output_file',
		help='The file where the optimized charges should be written to.',
		default='fitted_point_charges.csv', metavar='fitted_point_charges.csv')

	args = parser.parse_args()

	return args


def main():
	'''Read the constraints, transform them into matrix form, 
	and then use them to fit the point charges.'''
	# Read command line arguments
	args = parse_command_line()

	# Look up the relationship between ASE indices, atom names
	pmd_struct, pmd_top, ase2pmd = create_structure(args.pdb_infile, args.top_infile)

	# Calculate constraints
	logic_constraints, charge_constraints = get_constraints(args, ase2pmd=ase2pmd)

	# Import A and B matrices from HORTON
	A, B = read_horton_cost_function(args.horton_cost_function)

	# Run the constrained minimization
	q, f = constrained_minimize(A, B, logic_constraints, charge_constraints)

	print(q)

	# Save charges
	write_charges(q)

	# Save Lagrange forces
	write_forces(f)

	print('Done.')
	

if __name__ == '__main__':
	main()
