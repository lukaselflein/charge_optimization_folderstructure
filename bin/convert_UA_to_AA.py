""" Change structure with implicit Hydrogen to one with explicitely defined H-atoms.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Johannes Hoermann <johannes.hoermann@imtek.uni-freiburg.de>
Modified: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import ase.io
import sys
import warnings
import numpy as np
import parmed as pmd

from ase.data import atomic_numbers
from ase.neighborlist import NeighborList
from matscipy.neighbours import neighbour_list
from parmed import gromacs
from smamp.insertHbyList import insertHbyList

def read_input_files(input_dir='../0_initial_structure'):
    """
    Search for and read input files (with implicit H-atoms).
    """
    ase_struct, pmd_top = None, None

    files = os.listdir(input_dir)
    for filename in files:
        if ('snapshot' in filename) and ('.pdb' in filename):
            ase_struct = ase.io.read(os.path.join(input_dir, filename))
            pmd_struct = pmd.load_file(os.path.join(input_dir, filename))
        elif '.top' in filename:
            pmd_top = gromacs.GromacsTopologyFile(os.path.join(input_dir, filename), 
                                                      parametrize=False)

    # Make sure we actually found everything we need
    if ase_struct is None:
        raise RuntimeError('structure file (.pdb) not found in {}'.format(input_dir))
    if pmd_top is None:
        raise RuntimeError('topology file (.top) not found in {}'.format(input_dir))
    return ase_struct, pmd_struct, pmd_top


def main():
    """Execute everything."""
    print('This is {}.'.format(__file__))
    # Read the united-atoms files extracted from the MD-simulation trajectory
    # throws some warnings on angle types, does not matter for bonding info
    with warnings.catch_warnings():
         warnings.simplefilter('ignore')
         ase_struct, pmd_struct, pmd_top = read_input_files()

    # define a dictionary, marking how many H atoms are missing at which bonding partner explicitly
    implicitHbondingPartners={'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}
    
    pmd_top.strip(':SOL,CL') # strip water and electrolyte from system

    pmd_top.box = pmd_struct.box # Needed because .prmtop contains box info
    pmd_top.positions = pmd_struct.positions

    # Insert the explicit hydrogens
    print('Inserting explicit hydrogens, please wait ...')
    with open('insert_H.log', 'w') as logfile:
       new_ase_struct, new_pmd_top = insertHbyList(ase_struct,pmd_top,implicitHbondingPartners,1.0, logfile)

    # Write output
    new_ase_struct.write('ase_pdbH.pdb')
    new_ase_struct.write('ase_pdbH.traj')

    # Write other output
    new_pmd_top.write_pdb('pmd_pdbH.pdb')
    test_pmd = pmd.load_file('pmd_pdbH.pdb')
    # some topology format, un functionality similar to GROMACS' .top, but readable by VMD
    # new_pmd_top.write_psf('pmd_pdbH.psf')
    print('Done.')


if __name__ == '__main__':
    main()
