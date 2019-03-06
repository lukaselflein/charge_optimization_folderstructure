# Initialize the top level folder structure
cp -r .template_simulation 1_charge_cycle

# Copy the GROMACS simulation output
# Needed: .top, .tpr, .xtc
mv simulation_data/example.top 1_charge_cycle/simulation
mv simulation_data/example.tpr 1_charge_cycle/simulation
mv simulation_data/example.xtc 1_chareg_cycle/simulation

# Create the snapshots and subfolder structure
cd 1st_run
module purge
module load gromacs/2016.4-gnu-5.2
python create_snapshots_from_trajectory.py -tpr simulation/example.tpr -top simulation/example.top -xtc simulation/example.xtc -s 100 -d 100 -e 1000

# Convert UA to AA
module load matscipy/0.2.0
python loop_convert_UA_to_AA.py
# Output: ase_pdbH.traj, (ase_pdbH.pdb, pmd_pdbH.pdb)


# Submit gpaw optimization, ESP & Rho job
module purge
module load devel/python/3.6.5
module load smamp
python loop_submit.sh
# Output: esp.cube, rho.cube, (rho_pseudo.cube)

# Convert AA to UA
module purge
module load gromacs
module load matscipy/0.2.0
module load smamp
python loop_convert_UA_to_AA.py
# Output: esp_ua.cube, rho_ua.cube

# Calculate cost functions with HORTON
module purge
module load horton/2.1.0b3
module load smamp
loop_cost_functions.py
# Output: cost.h5

# Average over cost functions
average_cost.py 

# Fit ESP cost function
# Output: test.top, fitted_point_charges.top
python fitESPconstrained.py ua_cost.h5 snapshot1.pdb example.top atoms_in_charge_group.csv charge_group_total_charge.csv atoms_of_same_charge.csv test.top fitted_point_charges.csv

# Write new point charges into .rtf topology file
#TODO

# Plot point charges
#TODO
