#!/bin/bash

# Post-processing steps to extract charges from finished DFT calculations

# Convert AA to UA
module purge
module load gromacs
module load matscipy/0.2.0
module load smamp
python loop_convert_AA_to_UA.py
# Output: esp_ua.cube, rho_ua.cube

# Calculate cost functions with HORTON
module purge
module load horton/2.1.0b3
module load smamp
python loop_cost_functions.py
# Output: cost.h5

# Average over cost functions
module purge
module load horton/2.1.0b3
module load smamp
python average_cost.py 
# Output: costfunction_average.h5

# Fit ESP cost function
module purge
module load horton/2.1.0b3
module load smamp
module load gromacs
python ../bin/fitESPconstrained.py horton_charges/costfunction_average.h5 500_ps_snapshot/0_initial_structure/snapshot.pdb 500_ps_snapshot/0_initial_structure/example.top ../fitting_constraint_files/atoms_in_charge_group.csv ../fitting_constraint_files/charge_group_total_charge.csv ../fitting_constraint_files/atoms_of_same_charge.csv updated_charges.top horton_charges/fitted_point_charges.csv --qtot 0 --verbose 
# Output: test.top, fitted_point_charges.top

# Write new point charges into .rtp topology file
module purge
module load devel/python/3.6.5
module load smamp
python charges_to_rtp.py -rtp md_simulation/n7nh2.rtp -csv horton_charges/fitted_point_charges.csv -out modified_charges_topology.rtp
# Output: modified_charges_topology.rtp

# Calculate indiviudal charges for error bars
module purge
module load devel/python/3.6.5
module load gromacs
module load smamp
python loop_charges.py
# Output: fitted_point_charges.csv

# Plot charges
module purge
module load devel/python/3.6.5
python plot_charges.py
# Output: scatter_q_unconstrained.png, scatter_before_after.png (in pictures/ )
