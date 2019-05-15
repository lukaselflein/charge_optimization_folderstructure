""" Test functions of the refactored fitESPconstrained.py module. 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import unittest
import numpy as np
import pandas as pd

from refactored_fitESPconstrained import unconstrained_minimize
from refactored_fitESPconstrained import constrained_minimize
from refactored_fitESPconstrained import symmetry_groups_to_matrix
from refactored_fitESPconstrained import make_symmetry_constraints
from refactored_fitESPconstrained import make_group_constraints
from refactored_fitESPconstrained import make_atom_name_constraints


class TestUnconstrainedMinimize(unittest.TestCase):
   """Testing unconstrained_minimize"""
   
   def test_1_d(self):
      """ Quadratic expression:
      0.5 * x^2 - x
      Find maximum, derivative is zero:
      x - 1 = 0
      x = 1
      """
      A = [[1]]
      B = [1]
      x = unconstrained_minimize(A=A, B=B)
      self.assertAlmostEqual(np.dot(A, x), B)

   def test_2_d(self):
      """
      3 * x0 + x1 = 9
      x0 + 2 * x1 = 8
      """
      A = np.array([[3, 1], [1, 2]])
      B = [9, 8]
      x = unconstrained_minimize(A=A, B=B)
      np.testing.assert_allclose(np.dot(A, x), B)


class TestConstrainedMinimize(unittest.TestCase):
   """Testing constrained_minimize"""

   def test_symmetry_2D(self):
      """
      3 * x0 + x1 = 9
      x0 + 2 * x1 = 8
      Symmetry Constraint: 
      x0 = x1
      """
      A = np.array([[3, 1], [1, 2]])
      B = np.array([9, 8])
      # Constraint: x0 and x1 are equal
      D = np.array([[1, -1]])
      Q = np.array([0])
      x, f = constrained_minimize(A=A, B=B, D=D, Q=Q)

      # Is the constraint fulfilled?
      self.assertEqual(x[0], x[1])

      # Is x an actual solution to the problem?
      A_con = np.array([[3, 1, 1], [1, 2, -1], [1, -1, 0]])
      B_con = np.array([9, 8, 0])
      x_lag = np.concatenate((x, f))
      np.testing.assert_allclose(np.dot(A_con, x_lag), B_con)


class TestSymmetryToMatrix(unittest.TestCase):
   """Testing symmetry_groups_to_matrix"""
   def test_small(self):
      groups = [[0, 2, 3]]
      attempt = symmetry_groups_to_matrix(groups, n_atoms=5)[0]
      correct_matrix = np.array([[ 1,  0, -1,  0,  0],
                              [ 1,  0,  0, -1,  0]])
      np.testing.assert_allclose(attempt, correct_matrix)

   def test_multiple_residue(self):
      pass


class TestGroupConstraints(unittest.TestCase):
   """Tests for make_group_constraints"""

   def test_simple_case(self):
      charge_groups = {1: [0, 1, 2], 2: [3, 4, 5]}
      group_q = pd.DataFrame([2, 5], index=[1, 2])
      n_atoms = 7
      D_matrix, Q_vector = make_group_constraints(charge_groups, group_q, n_atoms)
      
      # The charge constraints should be in the Q_vector
      np.testing.assert_allclose(Q_vector, [2, 5])

      # The logic constraints should be of the form [1, 1... 0]
      correct_constraints = np.array([[1, 1, 1, 0, 0, 0, 0],
                                      [0, 0, 0, 1, 1, 1, 0]])
      np.testing.assert_allclose(D_matrix, correct_constraints)

   def test_multiple_residue(self):
      """Are charge groups defined for each residuum seperately?"""

   def test_different_charges(self):
      """Do groups get assigned their respective charges?"""
      pass


class TestAtomNameConstraints(unittest.TestCase):
   """Tests for make_atom_name_constraints(ase2pmd)"""

   def simple_dict(self):
      """Are atoms of same name constrained to have equal charge?"""
      ase2pmd = {0: ('CA1', 'terA'), 1: ('CA1', 'terB'), 2: ('CA1', 'terC')}
      correct_constraints = np.array([[ 1, -1,  0],
                                      [ 1,  0, -1]]) 
      D_matrix, Q_vector = make_atom_name_constraints(ase2pmd)
      np.testing.assert_allclose(D_matrix, correct_constraints)


if __name__ == '__main__':
    unittest.main()
