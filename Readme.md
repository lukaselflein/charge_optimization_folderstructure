# Sarah Folder Structure
Create a consistent, documented folder structure for charge optimization.

### Organization
The folders are organized in four levels:
1. Lowerest: subfolders for the diverse charge fitting and conversion routines.
2. Snapshot: Each snapshot of has its own folder
3. Simulation: For each charge, a seperate MD simulation is run which is divided into snapshots
4. Outer loop: Simulating and fitting charges is iteratated until charges converge. This is the highest level.
