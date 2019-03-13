""" Extract atomic charges from Bader Charge Analysis output.
Copyright 2019 Simulation Lab
University of Freiburg
Author:
Modified: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import ase
from parmed import gromacs
import numpy as np
from smamp.insertHbyList import insertHbyList
import warnings

def main():
	""" Execute everything."""
	# Read the atomic structure
	ase_struct = ase.io.read('snapshot.pdb')

	# Read the toplogy
	with warnings.catch_warnings():  # supress warnings
		warnings.simplefilter("ignore")
		pmd_top = gromacs.GromacsTopologyFile('example.top',parametrize=False)

	implicitHbondingPartners={'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}
	new_ase_struct, new_pmd_top, names, residues = insertHbyList(ase_struct, pmd_top, 
								     implicitHbondingPartners)

	atomic_number = new_ase_struct.get_atomic_numbers()
	bader_charges = np.genfromtxt('ACF.dat',skip_header=2,skip_footer=4,usecols=4)
	atomic_charges = atomic_number - bader_charges  # this inverts the sign of the bader charges

	results =[]
	for i in range(127):
		results.append((names[i], residues[i], atomic_charges[i],
				atomic_number[i],bader_charges[i]))

	# find the summed charges of carbon atoms and corresponding H atoms
	list_of_C_atoms = ['CD3', 'CD4', 'CA2', 'CA3', 'CB2', 'CB3']
	atom_res_totq = []  # list of [C-atom, residue, total_charge]
	nra = np.array([names,residues,atomic_charges])

	total_charge_results = np.array(results)[:,:4]
	total_charge_results[:,3] = total_charge_results[:,2]

	terminal_names = np.unique(nra[1])
	for term in terminal_names:
		residue = np.where(nra[1] == str(term))
		nra_for_one_residue = nra[:,residue[0]]

		# find all carbon atoms and the belonging H atoms
		for c_atom in list_of_C_atoms:
			C_H_positions = np.flatnonzero(np.core.defchararray.find(
							nra_for_one_residue[0], str(c_atom))!=-1)
			C_H_atoms = nra_for_one_residue[:, C_H_positions]
			total_charge = (C_H_atoms[2]).astype(np.float).sum()
			atom_res_totq.append([c_atom, term, total_charge])

			#overwrite the charge of list_of_C_atoms in results with total_charge
			res_position = np.where((total_charge_results[:,0]==c_atom) 
						& (total_charge_results[:,1]==term))
			#print(res_position+(2,))
			total_charge_results[res_position+(3,)] = total_charge

	with open('bader_charges.csv','w') as f:
		f.write("atom,\tresidue,\tbader charge"+ '\n')
		for line in total_charge_results:
			f.write(', \t '.join(line[:-1])+'\n')

if __name__ == '__main__':
	main()
