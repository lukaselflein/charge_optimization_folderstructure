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
# Output: ase_pdbH.traj, ase_pdbH.pdb, pmd_pdbH.pdb

# Submit gpaw optimization, ESP & Rho job
python loop_submit_gpaw.sh
# Output: esp.cube, rho.cube

# Convert AA to UA
module purge
module load numlib/python_scipy/0.19.2-python_numpy-1.11.2-python-2.7.12
module load gromacs
module load smamp
python loop_convert_AA_to_UA.py
# python aa2ua_cube.py ../0_initial_structure/snapshot.pdb  ../0_initial_structure/example.top  ../2_dft_calculations/esp.cube esp_ua.cube 2>&1 > esp_conversion.log
# python aa2ua_cube.py ../0_initial_structure/snapshot.pdb  ../0_initial_structure/example.top  ../2_dft_calculations/rho.cube rho_ua.cube  2>&1 > rho_conversion.log
# Output: esp_ua.cube, rho_ua.cube

