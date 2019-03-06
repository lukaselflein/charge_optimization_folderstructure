# First, we need gromacs for splitting the trajectory
module purge
module load gromacs/2016.4-gnu-5.2

# Create the snapshots and subfolder structure
python create_snapshots_from_trajectory.py -tpr simulation/example.tpr -top simulation/example.top -xtc simulation/example.xtc -s 100 -d 100 -e 1000

# We need scipy for the next step
module load devel/python/3.6.5
module load smamp

# Convert UA to AA
python loop_convert_UA_to_AA.py 
# Output: ase_pdbH.traj, (ase_pdbH.pdb, pmd_pdbH.pdb)

# Submit gpaw optimization, calc Electrostatic potential and density
python loop_submit_gpaw.sh
# Output: esp.cube, rho.cube, (rho_pseudo.cube)
