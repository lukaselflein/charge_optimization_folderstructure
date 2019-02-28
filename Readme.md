# Sarah Folder Structure
Create a consistent, documented pipline and folder structure for self-consistent MD charge optimization.

### Prerequisites
Make sure you have the smamp module:
1. Local installation
pip install --user smamp

2. Module installation
Alternatively, you can use the make_module.sh script to install & load it as a module.

3. Use version from Lukas
module use /home/fr/fr_fr/fr_le1019/modulefiles
module load smamp

### Organization
The folders are organized in four levels:
1. Lowerest: subfolders for the diverse charge fitting and conversion routines.
2. Snapshot: Each snapshot has its own folder, containing (1).
3. Simulation: For each set of charges, a seperate MD simulation is run. Contains subfolders for the snapshots extracted from the trajectory (2).
4. Current directory: Simulating and fitting charges is iteratated until charges converge. This is the highest level. Subfolders correspond to different simulations (3).

### Content
`bin`: Bash and python scripts for fitting, conversion, extracting, and visualization
`simulation_template`: The full template simulation folder structure.

### Usage
Copy the `simulation_template` folder structure, and rename to a specific simulation name, e.g. `simulation-1.0-a`. 

Then, run the simulations, make sure the data is saved in the appropriate location.
