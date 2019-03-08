#!/bin/bash

# Preprocessing and submission steps for charge optimization

# Initialize the top level folder structure
pwd
cp -r .template_simulation 1_charge_cycle
cd 1_charge_cycle
pwd

# Copy the GROMACS md_simulation output
# Needed: .top, .tpr, .xtc
mkdir md_simulation
cp ../md_simulation_data/example.top md_simulation/
cp ../md_simulation_data/example.tpr md_simulation/
cp ../md_simulation_data/example.xtc md_simulation/
cp ../md_simulation_data/n7nh2.rtp   md_simulation/

# Create the snapshots and subfolder structure
module purge
module load gromacs/2016.4-gnu-5.2
python create_snapshots_from_trajectory.py -tpr md_simulation/example.tpr -top md_simulation/example.top -xtc md_simulation/example.xtc -s 100 -d 100 -e 1000

# Convert UA to AA
module purge
module load gromacs/2016.4-gnu-5.2
module load matscipy/0.2.0
python loop_convert_UA_to_AA.py
# Output: ase_pdbH.traj, (ase_pdbH.pdb, pmd_pdbH.pdb)


# Submit gpaw optimization, ESP & Rho job
module purge
module load devel/python/3.6.5
module load smamp
python loop_submit.py
# Output: esp.cube, rho.cube, (rho_pseudo.cube)

