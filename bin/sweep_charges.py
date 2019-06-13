""" Calculate HORTON charges for all .cube files in folderstructure.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import sys
import shutil
import smamp
import warnings
import subprocess

from smamp.tools import cd
from smamp.tools import read_total_charge
from smamp.tools import find
from fitESPconstrained import create_structure, read_horton_cost_function
from fitESPconstrained import get_constraints, constrained_minimize, unconstrained_minimize
from fitESPconstrained import write_charges, write_forces


def calc_charges(pdb_infile, top_infile, hydrogen_file, horton_cost_function, output_file):
   '''Wraps fitESPconstrained.py'''
   # Look up the relationship between ASE indices, atom names
   pmd_struct, pmd_top, ase2pmd = create_structure(pdb_infile, top_infile,
                                                   hydrogen_file)
   print('Atomic structure built.')

   # Import A and B matrices from HORTON
   A, B = read_horton_cost_function(horton_cost_function)

   charge_group_file = None 
   charge_group_charges_file = None
   symmetry_file = None

   # Calculate constraints
   logic_constraints, charge_constraints = get_constraints(args=None, ase2pmd=ase2pmd, 
                                                    charge_group_file=charge_group_file, 
                                                    charge_group_charges_file=charge_group_charges_file,
                                                    symmetry_file=symmetry_file,
                                                    debug=False)

   print('Constraints caluclated: {} non-redunant.'.format(logic_constraints.shape[0]))

   # Run the constrained minimization
   q, f = constrained_minimize(A, B, logic_constraints, charge_constraints)
   # print('Constrained minimization done.')
   # print('Extremal charges: {:1.5f}, {:1.5f}'.format(q.min(), q.max()))
   # print('Extremal Lagrange forces: {:1.5f}, {:1.5f}'.format(f.min(), f.max()))

   q_unconstrained = unconstrained_minimize(A, B)

   # Save charges
   charge_df = write_charges(q, q_unconstrained, ase2pmd, out_name=output_file, plot=False)

   # Save Lagrange forces
   # write_forces(f, logic_constraints, ase2pmd)
   print('Charges and forces written.')

def main():
	""" Execute everything."""
	print('This is {}.'.format(__file__))

	# We need any one top and pdb file for the ordering of atom-names;
	# The exact snapshot we use does not matter as the ordering is invariant
        pdb_file = find(path='.', folder_keyword='0_initial_structure', file_keyword='.pdb')[0]
        top_file = find(path='.', folder_keyword='0_initial_structure', file_keyword='.top')[0]
        hyd_file = find(path='..', folder_keyword='fitting', file_keyword='hydrogen_per_atom.csv')[0]

	# Also, we need constraint files

	for lnrho in range(-10, 0):
		print('lnrho = {} .'.format(lnrho))
		# Find the path for the average cost function
		cost_avg = find(path='.', folder_keyword='horton_charges/sweep_rhoref', 
		    	   file_keyword='costfunction_average_{}.h5'.format(lnrho), 
	 		   nr_occ=1)[0]

		folder = os.path.split(cost_avg)[0]
		output_file = os.path.join(folder, 'avg_charges_{}.csv'.format(lnrho))

		calc_charges(pdb_file, top_file, hyd_file, cost_avg, output_file=output_file)
	print('Done.')
	exit()

	for lnrho in range(-10, 0):

		# Find the paths for the unaveraged snapshot cost functions
		cost_paths = find(path='.', folder_keyword='4_horton_cost_function/lnrho_sweep', 
		    	          file_keyword='cost_{}.h5'.format(lnrho), 
     	 		          nr_occ=None)

		print(cost_paths)


if __name__ == '__main__':
	main()
