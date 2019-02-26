#!/usr/bin/env python
""" Minimize the enegy of input structure and extract ESP. 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

from ase.io import read
from ase.io import write

from gpaw import GPAW
from gpaw import restart

from gpaw import GPAW
from ase.optimize.bfgslinesearch import BFGSLineSearch #Quasi Newton
from ase.units import Bohr
from ase.units import Hartree

import os.path
import argparse
import io
import numpy as np


def parser():
	"""
	Parse Command line arguments.
	"""
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-t', '--trajectory', metavar='ase_pdbH.traj', default='ase_pdbH.traj', help='The path to the trajectory file.')
	parser.add_argument('-r', '--restart', metavar='restart.gpw', default=None, help='The path to a restart file, optional.')

	args = parser.parse_args()
	traj_file = args.trajectory
	gpw_file = args.restart

	return traj_file, gpw_file

def minimize_energy(traj_file):
	"""
	Run a BFGS energy minimization of the smamp molecule in vacuum.

	Args:
	traj_file: the path to a trajectory file (all atom format)

	Returns:
	struc: an internal molecular structure object
	calc: internal calculation object	
	"""

	# Read in the trajectory file
	struc = read(traj_file)
	# Set up the box
	struc.set_cell([25,25,25])
	struc.set_pbc([0,0,0])
	struc.center()
	# Define gpaw convergence&simulation parameters
	calc  = GPAW(xc='PBE', h=0.2, charge=0,
		     spinpol=True, convergence={'energy': 0.001})
	struc.set_calculator(calc)
	dyn = BFGSLineSearch(struc, trajectory='molecule.traj',
			     restart='bfgs_ls.pckl', logfile='BFGSLinSearch.log')

	# run the simulation
	dyn.run(fmax=0.05)

	# Maybe this is useful? Does not seem to be needed.
	# Epot  = struc.get_potential_energy()

	# Save everything into a restart file
	calc.write('restart.gpw', mode='all')

	return struc, calc

def read_restart(gpw_file):
	""" Extract the structure and calculation from a restart file."""
	struc, calc = restart(gpw_file)
	return struc, calc

def extract(struc, calc):
	"""
	Extract & write electrostatic potential and densities.
	
	Arg:
	struc: an internal molecular structure object
	calc: internal calculation object
	"""
	# Extract the ESP
	esp = calc.get_electrostatic_potential()

	# Convert units
	esp_hartree = esp / Hartree   
	write('esp.cube', struc, data=esp_hartree)

	# Psedo-density, what does this do?
	rho_pseudo      = calc.get_pseudo_density()
	rho_pseudo_per_bohr_cube = rho_pseudo * Bohr**3
	write('rho_pseudo.cube', struc, data=rho_pseudo_per_bohr_cube) 

	# Density
	rho             = calc.get_all_electron_density()
	rho_per_bohr_cube = rho * Bohr**3
	write('rho.cube', struc, data=rho_per_bohr_cube) 

def main():
	"""
	Execute everything.
	"""

	# Read Command line arguments	
	traj_file, gpw_file = parser()

	# Check if a restart file was provided
	if gpw_file is not None:
		# If we have a restart file, use everything from it
		struc, calc = read_restart(gpw_file)
	
	# Otherwise, we need to optimize based on our input file first.
	else:
		struc, calc = minimize_energy(traj_file)

	# Now we can extract ESP, Rho, etc.
	extract(struc, calc)
	print('Done.')


if __name__ == '__main__':
	main()
