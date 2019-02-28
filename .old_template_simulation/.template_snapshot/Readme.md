# Template for the charge analysis of one snapshot

This is the folder structure to use in doing charge optimization.
In each of these folders, one snapshot of a MD-simulation is contained, the relevant molecule extracted, and best-fit charges are determined. These are then written into the input file for a new MD simulation.

### Subfolders:
* `0_initial_structure`: the MD snapshot, extraced topology and xyz files (UnitedAtoms)
* `1_all_atoms_structure`: the inital structure converted to the AllAtoms format
* `2_dft_calculations`: gpaw files determining the electronic (DFT) structure of the molecule
* `3_electron_densities`: electron density and ESP extracted from the gpaw calculations
* `4_horton_charge_fitting`: charges determined by HORTON (fitting an ESP cost-function)
* `5_bader_charge_fitting`: charges determined via Bader (charges within zero flux surfaces)

### Flow of information
Information flows in the direction of the subfolders as listed above.
A UnitedAtoms structure is converted to AllAtoms, the corresponding density calculated and optimal charges are determined.
In the folder structure above, the charges are then compared for the different snapshots.
In the hierarchy level yet one above, the average charges are used to run another MD simulation, and the time-evolution of the charges is analyzed - hopefully, they converge to a stationry value, which can then be used for production-level MD.
