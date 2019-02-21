#!/usr/bin/env python
""" Search for inital united atom structures, convert them to all atom format.
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
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-t', '--trajectory', metavar='ase_pdbH.traj', default='ase_pdbH.traj')
	parser.add_argument('-r', '--restart', metavar='ase_pdbH.traj', default='ase_pdbH.traj')

	args = parser.parse_args()
	traj_file = args.infile_traj
	gpw_file  = args.outfile_gpw
	box       = args.box

	return traj_file, gpw_file, charge, box

def minimize_energy(traj_file):
	"""
	Run a BFGS energy minimization of the smamp molecule in vacuum.
	"""
	struc = read(traj_file)
	struc.set_cell([25,25,25])
	struc.set_pbc([0,0,0])
	struc.center()
	calc  = GPAW(xc='PBE', h=0.2, charge=0,
		     spinpol=True, convergence={'energy': 0.001})

	struc.set_calculator(calc)
	dyn = BFGSLineSearch(struc, trajectory='molecule.traj',
			     restart='bfgs_ls.pckl', logfile='BFGSLinSearch.log')
	dyn.run(fmax=0.05)

	# Maybe this is useful? Does not seem to be needed.
	Epot  = struc.get_potential_energy()

	# Save everything into a restart file
	calc.write('restart.gpw', mode='all')

	return struc, calc

def read_restart(gpw_file):
	struc, calc = restart(gpw_file)
	return struc, calc

def extract(struc, calc):
	"""
	Extracts & writes electrostatic potential and densities.
	"""
	# struc, calc = restart(gpw_file)

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

	
	traj_file, gpw_file = parser()
	if gpw_file is not None:
		struc, calc = read_restart
	
	else:
		traj_file = 'ase_pdbH.traj'
		struc, calc = minimize_energy(traj_file)

	extract(struc, calc)


if __name__ == '__main__':
	main()
