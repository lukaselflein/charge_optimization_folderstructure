""" Fit point charges to a HORTON costfunction under constraints.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
Based on legacy code by Johannes Hormann
"""

import argparse


def parse_command_line():
	"""Read file locations from command line interface."""
	parser = argparse.ArgumentParser(prog='esp-fit-constrained.py',
        		description='Estimate charges from a HORTON ESP cost function under constraints.')
	
        parser.add_argument('-hor', '--horton_cost_function',
		help='The location of the HORTON cost function file.',
		required=True, metavar='cost.h5')

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


def constrained_minimize(A, B, logic_constraints, charge_constraints):
	x = [0]
	f = 0
	return x, f


def parse_charge_groups(path):
	print(path)


def parse_group_charges(path):
	print(path)


def parse_symmetry(path):
	print(path)


def make_constraints(*args, **kwargs):
	pass


def remove_redunancy(*args, **kwargs):
	pass


def get_constraints(args):
	'''Read provided constraint files and convert them into matrix form.'''
	charge_group_file = args.charge_groups
	charge_group_charges_file = args.charge_group_charges
	symmetry_file = args.symmetry_file
	
	if charge_group_file is not None:
		if charge_group_charges_file is None:
			err = 'Charge groups defined: {}'.format(charge_group_file)
			err += '\n But no total charges were defined.'
			raise ValueError(err)

		charge_groups = parse_charge_groups(charge_group_file)
		charges = parse_group_charges(charge_group_charges_file)
		group_constraints = make_constraints(charge_groups, charges)

	else:
		group_constraints = None

	if symmetry_file is not None:
		symmetry = parse_symmetry(symmetry_file)
		symmetry_constraints = make_constraints(symmetry_file)
	else:
		symmetry_constraints = None

	# Remove redunant constraints
	reduced_constraints = remove_redunancy(group_constraints,
					       symmetry_constraints)
	
	logic_constraints, charge_constraints = None, None
	return logic_constraints, charge_constraints


def import_horton(path):
	print(path)
	A = None
	B = None
	return A, B


def write_charges(charges):
	pass


def write_forces(f):
	pass


def get_ase_indices():
	pass

def main():
	'''Read the constraints, transform them into matrix form, 
	and then use them to fit the point charges.'''
	# Read command line arguments
	args = parse_command_line()

	# Look up the relationship between ASE indices, atom names
	get_ase_indices()

	# Calculate constraints
	logic_constraints, charge_constraints = get_constraints(args)

	# Import A and B matrices from HORTON
	A, B = import_horton(args.horton_cost_function)

	# Run the constrained minimization
	q, f = constrained_minimize(A, B, logic_constraints, charge_constraints)

	# Save charges
	write_charges(q)

	# Save Lagrange forces
	write_forces(f)

	print('Done.')
	

if __name__ == '__main__':
	main()
