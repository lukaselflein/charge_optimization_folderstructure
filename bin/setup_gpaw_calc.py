#setup the gpaw calculation

from ase.io import read
from gpaw import GPAW
#from ase.optimize import FIRE #Quasi Newton + friction
from ase.optimize.bfgslinesearch import BFGSLineSearch #Quasi Newton

traj_file = '../1_all_atom_structure/ase_pdbH.traj'
print('Reading {}'.format(traj_file))
struc = read(traj_file)
struc.set_cell([25,25,25])
struc.set_pbc([0,0,0])
struc.center()
calc  = GPAW(xc='PBE', h=0.2, charge=0,
             spinpol=True, convergence={'energy': 0.001})

struc.set_calculator(calc)
#opt   = FIRE(struc, trajectory='molecule.traj', logfile='fire.log')
dyn = BFGSLineSearch(struc, trajectory='molecule.traj',
                     restart='bfgs_ls.pckl', logfile='BFGSLinSearch.log')
dyn.run(fmax=0.05)

# We can actually just write out a restart file
gpw_file = 'system.gpw'
calc.write(gpw_file, mode='all')

