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
import sympy

import parmed as pmd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from smamp.insertHbyList import insertHbyList
from smamp.tools import read_atom_numbers


def create_structure(infile_pdb, infile_top, hydrogen_file, strip_string=':SOL,CL'):
    """Build ase-format atomic structure descriptions.
    Especially useful is a dictionary listing the relationship between ase indices and atom names.

    Args:
       infile_pdb (str): path to the gromacs structure file
       infile_top (str): path to the gromacs topology file
       hydrogen_file (str): file with explicit hydrogen atom description
       strip_string (str): atoms to be removed from .pdb file

    Returns:
       pmd_struct:
       pmd_top: 
       ase2pmd (dict): A map of ase indices to atom names
    """
        
    implicitHbondingPartners = read_atom_numbers(hydrogen_file)
        
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

def constrained_minimize(A, B, D=None, Q=None):
   """Find the minimum of the HORTON cost function.
   The cost function is parametrized with matrix A and vector B.
   In the unconstrained case, the minimization is equivalent to solving
   A x - B = 0
   for the charges x.

   In the case of constraints, we have to solve the problem
   A D  x  =  B  
   D 0  l     Q
   with D being the logical constraints, and Q the respective charge values.

   The function first stacks (A, D) and (B, Q) to resemble 
   the unconstrained case formally, then solves the constrained equation.

   Args:
      A (np.array): Matrix with quadratic terms of cost fucnction
      B (np.array): Vector with linear tearms of cost function
      D (np.array): Matrix with constraint logic
      Q (np.array): Vector with constraint charges

   Returns:
      charges (np.array): Vector of optimal charges
      langrange_forces (np.array): Vector of forces neccesary constrain charges
   """
   # Default to zero total charge constraint
   if D is None and Q is None:
      Q = np.array([0])
      D = np.ones(B.shape[0])
      
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
   zeros = np.zeros((Q.shape[0], Q.shape[0]))
   A_con = stack([[A, D.T],
                  [D, zeros]])
   B_con = stack([B, Q]).T

   x = np.linalg.solve(A_con, B_con)

   charges = x[:len(B)]
   lagrange_forces = x[len(B):]
   return charges, lagrange_forces


def unconstrained_minimize(A, B):
   """Find the unconstrained minimum of the HORTON cost function A x - B = 0.

   Args:
      A (np.array): Matrix with quadratic terms of cost fucnction
      B (np.array): Vector with linear tearms of cost function

   Returns:
      charges (np.array): Vector of optimal charges
   """
   charges = np.linalg.solve(A, B)
   return(charges)


def parse_charge_groups(file_name, ase2pmd):
   """Read the charge group file."""
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
   """Read the file specifying total charges of each charge group."""
   group_q = pd.read_csv(file_name, sep=',', header=None, comment='#',
                         names=['charge'], index_col=0)
   group_q.charge = group_q.charge.astype(int)    
   return group_q

def parse_symmetry(file_name):
   """Read the file containing pair-symmetry constraints."""
   df = pd.read_csv(file_name, sep=',', header=None,  comment='#')
   symm_names = df.values.tolist()
   return symm_names
                
def symmetry_names_to_index_groups(symm_names, ase2pmd):
   """Transform atom-name based constraints to index-based constraints."""
   symm_groups = []
   for i in range(len(symm_names)):
      names = symm_names[i]
      symm_groups += [[]]
      for ase_index, atom_residuum in ase2pmd.items():
         # If the atom names match, pick the ase index 
         atom_name = atom_residuum[0]
         if names[0] == atom_name:
            # Every member of this group is supposed to have equal charge           
            symm_groups[i] += [ase_index]
         if names[1] == atom_name:
            symm_groups[i] += [ase_index]
   return symm_groups
            
def symmetry_groups_to_matrix(symm_groups, n_atoms):
   """Generate matrix-constraints from groups of same-charge indices.
   >>> groups = [[0, 2, 3]]
   >>> symmetry_groups_to_matrix(groups, n_atoms=5)[0]
   array([[ 1,  0, -1,  0,  0],
          [ 1,  0,  0, -1,  0]])
   """
   symm_list = []
   for group in symm_groups:
      for atom_index in group[1:]:
         matrix_row = np.zeros(n_atoms, dtype=int)
         matrix_row[group[0]] = 1
         matrix_row[atom_index] = -1
         symm_list += [matrix_row]

   symmetry_matrix = np.array(symm_list)
   symmetry_q = np.zeros(symmetry_matrix.shape[0], dtype=int)
   return symmetry_matrix, symmetry_q

def make_symmetry_constraints(symm_names, ase2pmd):
   """Transform atom-name symmetry constraints to ase-index matrix format."""
   symm_groups = symmetry_names_to_index_groups(symm_names, ase2pmd)
   n_atoms = len(ase2pmd)
   D_matrix, Q_vector = symmetry_groups_to_matrix(symm_groups, n_atoms)
   return D_matrix, Q_vector


def make_group_constraints(charge_groups, group_q, n_atoms):
   """Transform atom-name group charge group constraints to ase-index matrix form."""
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


def make_atom_name_constraints(ase2pmd):
   """Construct constraints for atoms of same name to have equal charge across residues."""
   # Extract unique atom names
   unique_names = []
   for ase_index, atom_residuum in ase2pmd.items():
      if atom_residuum[0] not in unique_names:
         unique_names += [atom_residuum[0]]

   name_groups = {}
   for name in unique_names:
      name_groups[name] = []

   # At which indices do atom names occur?
   for name in unique_names:
      for ase_index, atom_residuum in ase2pmd.items():
         if name in atom_residuum:
            name_groups[name] += [ase_index]

   # Keep name-groups with at least two members, don't need the rest 
   groups = []
   for name, index_list in name_groups.items():
      if len(index_list) > 1:
         groups += [index_list]

   # Transform the groups to matrix form
   groups = np.array(groups)
   D_matrix, Q_vector = symmetry_groups_to_matrix(groups, n_atoms=len(ase2pmd))

   return D_matrix, Q_vector


def nonsingular_concat(X, vector):
   """Appends vector to matrix X iff the resulting matrix is nonsingular.

   Args: 
      X (np.array): NxM Matrix to be appended to
      vector (np.array): Nx1 vector to be appended to X

   Returns:
      new_X (np.array): Nx(M+1) Matrix or None
   """ 
   # Cast vector to matrix
   vector = np.atleast_2d(vector)
   # Append vector as new row at bottom of matrix
   new_X = np.concatenate((X, vector), axis=0)

   # Check if matrix is still non-singular
   if new_X.shape[0] == np.linalg.matrix_rank(new_X):
      return new_X
   else:
      return None


def stack_constraints(X, Q_x, Y, Q_y):
   """Transform two constraint matrices/vector pairs into a single pair.

   Args:
      X (np.array): Constraint matrix to be appended to
      Y (np.array): Constraint matrix to be conatenated
      Q_x (np.array): The constraint charges corresponding to X
      Q_y (np.array): Constraint charges corresponding to Y
   """

   con_matrix = X.copy()
   con_q = Q_x.copy()

   if all([obj is not None for obj in (X, Y, Q_x, Q_y)]):
      for row in range(Y.shape[0]):
         new_matrix = nonsingular_concat(con_matrix, Y[row, :])
         if new_matrix is not None:
            con_matrix = new_matrix
            con_q = np.concatenate((con_q, np.atleast_1d(Q_y[row])))
         else:
            with open('dropped_constraints.txt', 'ab') as outfile:
               np.savetxt(outfile, Y[row, :], fmt='%d', newline=" ")
               outfile.write(b'\n')

   return con_matrix, con_q

   if X is None and (Y is not None and Q_y is not None):
      return Y, Q_y

   if (X is not None and Q_x is not None) and Y is None:
      return X, Q_x

   else:
      raise ValueError('Invalid mixture of empty and non-empty constraints')


def get_constraints(args, ase2pmd):
   '''Read provided constraint files and convert them into matrix form.'''
   charge_group_file = args.charge_groups
   charge_group_charges_file = args.charge_group_charges
   symmetry_file = args.symmetry_file

   name_matrix, name_q = make_atom_name_constraints(ase2pmd)
   
   if charge_group_file is not None:
      if charge_group_charges_file is None:
         err = 'Charge groups defined: {}'.format(charge_group_file)
         err += '\n But no total charges were defined.'
         raise ValueError(err)

      charge_groups = parse_charge_groups(charge_group_file, ase2pmd)
      group_q = parse_group_charges(charge_group_charges_file)
      n_atoms = len(ase2pmd)
      group_matrix, group_q = make_group_constraints(charge_groups, group_q, n_atoms)

   else:
      group_matrix, group_q = None, None

   if symmetry_file is not None:
      symmetry = parse_symmetry(symmetry_file)
      symmetry_matrix, symmetry_q = make_symmetry_constraints(symmetry, ase2pmd)
   else:
      symmetry_matrix, symmetry_q = None, None

   group_symm_matrix, group_symm_q = stack_constraints(group_matrix, group_q,
                                                       symmetry_matrix, symmetry_q)
   constraint_matrix, constraint_q = stack_constraints(group_symm_matrix, group_symm_q,
                                                       name_matrix, name_q)

   np.savetxt('symm_matrix.txt', symmetry_matrix, fmt='%d')
   np.savetxt('name_matrix.txt', name_matrix, fmt='%d')
   np.savetxt('group_matrix.txt', group_matrix, fmt='%d')
   np.savetxt('constraint_matrix.txt', constraint_matrix, fmt='%d')
   return constraint_matrix, constraint_q


def read_horton_cost_function(file_name):
   """Extract A and B HORTON cost function matrics from HDF5 binary."""
   cost_function = h5py.File(file_name)
   A = cost_function['cost']['A'][()]
   B = cost_function['cost']['B'][()]
   return A, B


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
   parser.add_argument('-hyd', '--hydrogen_file',
      help='The hydrogen insertion rules',
      default='hydrogen_per_atom.csv', metavar='hydrogen_per_atom.csv')
   return parser.parse_args()


def write_charges(q, q_unconstrained, ase2pmd, out_name='fitted_point_charges', plot=False):
   """Write array of charges into .csv output file."""
   def number_to_atom_name(i):
       return ase2pmd[i][0]
   def number_to_residuum(i):
       return ase2pmd[i][1]

   df = pd.DataFrame(q, columns=['q'])
   df['q_unconstrained'] = q_unconstrained
   df['indices'] = df.index
   df['atom'] = df.indices.apply(number_to_atom_name)
   df['residue'] = df.indices.apply(number_to_residuum)
   df = df.drop(['indices'], axis=1)
   df = df[['atom', 'residue', 'q', 'q_unconstrained']]
   df.to_csv(out_name)

   if plot:
      plt.plot(q, range(len(q)), lw=0, marker='o')
      plt.plot(q_unconstrained, range(len(q_unconstrained)), lw=0, marker='o')
      plt.show()

   return df


def write_forces(forces, logic_constraints, ase2pmd):
   """Write lagrange forces to .csv output file."""

   force_constraint = []
   if logic_constraints is not None:
      forces = np.atleast_2d(forces).T
      c = np.concatenate((forces, logic_constraints), axis=1)
      c = c[c[:,0].argsort()]

      # Sorted forces and constraints
      forces = c[:, 0]
      l = c[:, 1:]

      for i in range(len(l)):
         line = l[i]
         f = forces[i]
         constraint = np.nonzero(line)[0]
         readable_con = [f]
         for number in constraint:
            atom = ase2pmd[number][0] + '/' + ase2pmd[number][1]
            readable_con += [atom]
         force_constraint += [readable_con]

   with open('lagrange_forces.csv', 'w') as outfile:
      outfile.write('force, atom names\n')
      for entry in force_constraint:
         line = '{0:.3f}, '.format(entry[0])
         line += ' '.join(entry[1:])
         line += '\n'
         outfile.write(line)


def main():
   '''Read the constraints, transform them into matrix form, 
   and then use them to fit the point charges.'''
   # Read command line arguments
   args = parse_command_line()

   print('Extracting structure via ASE ...')
   # Look up the relationship between ASE indices, atom names
   pmd_struct, pmd_top, ase2pmd = create_structure(args.pdb_infile, args.top_infile, 
                                                   args.hydrogen_file)
   print('Atomic structure built.')

   # Import A and B matrices from HORTON
   A, B = read_horton_cost_function(args.horton_cost_function)

   # Calculate constraints
   logic_constraints, charge_constraints = get_constraints(args, ase2pmd=ase2pmd)
   print('Constraints caluclated, {} remaining.'.format(logic_constraints.shape[0]))

   # Run the constrained minimization
   q, f = constrained_minimize(A, B, logic_constraints, charge_constraints)
   print('Constrained minimization done.')

   q_unconstrained = unconstrained_minimize(A, B)

   # Save charges
   charge_df = write_charges(q, q_unconstrained, ase2pmd, out_name=args.output_file, plot=False)

   # Save Lagrange forces
   write_forces(f, logic_constraints, ase2pmd)
   print('Charges and forces written.')
   print('Done.')
   

if __name__ == '__main__':
   main()
