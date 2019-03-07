# Initialize the top level folder structure
cp -r .template_md_simulation 1_charge_cycle

# Copy the GROMACS md_simulation output
# Needed: .top, .tpr, .xtc
mv md_simulation_data/example.top 1_charge_cycle/simulation
mv md_simulation_data/example.tpr 1_charge_cycle/simulation
mv md_simulation_data/example.xtc 1_charge_cycle/simulation
mv md_simulation_data/n7nh2.rtp   1_charge_cycle/simulation

# Create the snapshots and subfolder structure
cd 1st_run
module purge
module load gromacs/2016.4-gnu-5.2
python create_snapshots_from_trajectory.py -tpr md_simulation/example.tpr -top simulation/example.top -xtc simulation/example.xtc -s 100 -d 100 -e 1000

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
python ../../bin/fitESPconstrained.py costfunction_average.h5 ../500_ps_snapshot/0_initial_structure/snapshot.pdb ../500_ps_snapshot/0_initial_structure/example.top atoms_in_charge_group.csv charge_group_total_charge.csv atoms_of_same_charge.csv test.top fitted_point_charges.csv --qtot 0 --verbose 
# Output: test.top, fitted_point_charges.top

# Write new point charges into .rtp topology file
module purge
module load smamp
module load gromacs
python charges_to_rtp.py -rtp md_simulation/n7nh2.rtp -csv horton_charges/fitted_point_charges.csv -out modified_charges_topology.rtp
# Output: modified_charges_topology.rtp

# Plot point charges
module purge
module load devel/python/3.6.5
python plot_charges
# Output: scatter_q_unconstrained.png, scatter_before_after.png (in pictures/ )
