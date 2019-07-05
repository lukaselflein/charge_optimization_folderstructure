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


def calc_charges(pdb_infile, top_infile, hydrogen_file, horton_cost_function, output_file,
                 charge_group_file=None, charge_group_charges_file=None, symmetry_file=None):
   '''Wraps fitESPconstrained.py'''
   # Look up the relationship between ASE indices, atom names
   pmd_struct, pmd_top, ase2pmd = create_structure(pdb_infile, top_infile,
                                                   hydrogen_file)
   print('Atomic structure built.')

   # Import A and B matrices from HORTON
   A, B = read_horton_cost_function(horton_cost_function)

   if charge_group_file is None:
      try:
         charge_group_file = find('..', 'fitting_constraint_files', 'atoms_in_charge_group.csv')[0]
      except:
              charge_group_file = None

   if charge_group_charges_file is None:
      try:
         charge_group_charges_file = find('..', 'fitting_constraint_files', 
                            'charge_group_total_charge.csv')[0]
      except:
         charge_group_charges_file = None
   
   if symmetry_file is None:
      try:
         symmetry_file = find('..', 'fitting_constraint_files', 'atoms_of_same_charge.csv')[0]
      except:
         symmetry_file = None
      

   # Calculate constraints
   logic_constraints, charge_constraints = get_constraints(args=None, ase2pmd=ase2pmd, 
                                                    charge_group_file=charge_group_file, 
                                                    charge_group_charges_file=charge_group_charges_file,
                                                    symmetry_file=symmetry_file,
                                                    debug=False)

   print('Constraints calculated: {} non-redunant.'.format(logic_constraints.shape[0]))

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
   print('Charges and written to {}.'.format(output_file))

def main():
   """ Execute everything."""
   print('This is {}.'.format(__file__))

   # We need any one top and pdb file for the ordering of atom-names;
   # The exact snapshot we use does not matter as the ordering is invariant
   pdb_file = find(path='.', folder_keyword='0_initial_structure', file_keyword='.pdb')[0]
   top_file = find(path='.', folder_keyword='0_initial_structure', file_keyword='.top')[0]
   hyd_file = find(path='..', folder_keyword='fitting', file_keyword='hydrogen_per_atom.csv')[0]

   cost_paths = find(path='.', folder_keyword='4_horton_cost_function/lnrho_sweep', 
           file_keyword='cost',
           nr_occ=None)
   lnrho_range = []
   sigma_range = []
   for charge_file in cost_paths:
      # Parse parameters from filename
      lnrho, sigma = charge_file[-15:-3].split('_')[-2:]
      lnrho_range += [lnrho]
      sigma_range += [sigma]

   lnrho_range = set(lnrho_range)
   sigma_range = set(sigma_range)

   for lnrho in lnrho_range:
      print('lnrho = {} ...'.format(lnrho))
      for sigma in sigma_range:

         # Find the paths for the unaveraged snapshot cost functions
         cost_paths = find(path='.', folder_keyword='4_horton_cost_function/lnrho_sweep',
                 file_keyword='cost_{}_{}.h5'.format(lnrho, sigma), 
                 nr_occ=None)

         # Find the path for the average cost function
         cost_avg = find(path='.', folder_keyword='horton_charges/sweep_rhoref', 
               file_keyword='costfunction_average_{}_{}.h5'.format(lnrho, sigma), 
               nr_occ=None)
         print(cost_avg)
         cost_paths += cost_avg

         for cost_file in cost_paths:
            folder = os.path.split(cost_file)[0]
            output_file = os.path.join(folder, 'charges_{}_{}.csv'.format(lnrho, sigma))

            if os.path.exists(output_file):
               # print('{} exists. Skipping ahead.'.format(output_file))
               continue

            print('Optimizing charges for {}.'.format(cost_file[:18]))
            calc_charges(pdb_file, top_file, hyd_file, cost_file, output_file=output_file)
      

   print('Done.')


if __name__ == '__main__':
   main()
