# Sarah Folder Structure
Create a consistent, documented folder structure for charge optimization.

### Organization
The folders are organized in four levels:
1. Lowerest: subfolders for the diverse charge fitting and conversion routines.
2. Snapshot: Each snapshot has its own folder, containing (1).
3. Simulation: For each set of charges, a seperate MD simulation is run. Contains subfolders for the snapshots extracted from the trajectory (2).
4. Current directory: Simulating and fitting charges is iteratated until charges converge. This is the highest level. Subfolders correspond to different simulations (3).
