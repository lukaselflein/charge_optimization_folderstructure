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
from gpaw import FermiDirac
from gpaw import Mixer

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
	parser.add_argument('-c', '--charge', metavar='2.0', default=None, help='The total charge of the molecule.')

	args = parser.parse_args()
	traj_file = args.trajectory
	gpw_file = args.restart
	charge = args.charge

	return traj_file, gpw_file, charge


def read_total_charge(path='../fitting_constraint_files/total_charge.csv'):
	"""Determine total charge of molecule by reading file.
	Arguments
	path: Path to the table.

	Returns:
	charge: a float
	"""
	# The default charge is zero, use this if no file is found
	if not os.path.isfile(path):
		print('WARNING: No charge file found at {}. Using default charge of 0.0'.format(path))
		return 0.0

	with open(path) as charge_file:
		for line in charge_file:
			charge = float(line)
	return charge


def minimize_energy(traj_file, charge):
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
	struc.set_cell([30,30,30])
	struc.set_pbc([0,0,0])
	struc.center()
	# Define gpaw convergence&simulation parameters
	calc  = GPAW(xc='PBE', h=0.2, charge=charge,
		     spinpol=True, convergence={'energy': 0.001},
		     mixer=Mixer(beta=0.25, nmaxold=10, weight=1.0),
		     occupations=FermiDirac(width=0.1))
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

	# Psedo-density, does not seem to be used in workflow
	# rho_pseudo      = calc.get_pseudo_density()
	# rho_pseudo_per_bohr_cube = rho_pseudo * Bohr**3
	# write('rho_pseudo.cube', struc, data=rho_pseudo_per_bohr_cube) 

	# Density
	rho             = calc.get_all_electron_density()
	rho_per_bohr_cube = rho * Bohr**3
	write('rho.cube', struc, data=rho_per_bohr_cube) 


def main():
	"""Execute everything."""
	print('This is {}.'.format(__file__))
	# Read Command line arguments	
	traj_file, gpw_file, charge = parser()

	# Use provided charge, or fallback to loading the charge from file
	if charge is None:
		print('No charge provided in command line arguments. Reading from default file ...')
		charge = read_total_charge(path='../../../fitting_constraint_files/total_charge.csv')
	print('A total charge of {} is used.'.format(charge))

	# Check if a restart file was provided
	if gpw_file is not None:
		print('Restart file {} provided. Reading it ...'.format(gpw_file))
		# If we have a restart file, use everything from it
		struc, calc = read_restart(gpw_file)
	
	# Otherwise, we need to optimize based on our input file first.
	else:
		print('No restart file provided. Starting a new minimization ...')
		struc, calc = minimize_energy(traj_file, charge)

	print('Minimization finished. Extracting ESP and Density ...')
	# Now we can extract ESP, Rho, etc.
	extract(struc, calc)
	print('Done.')


if __name__ == '__main__':
	main()
