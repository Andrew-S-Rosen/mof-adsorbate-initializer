import os
from mai.adsorbate_constructor import adsorbate_constructor
from ase.io import read

grid_path = os.path.join('example_MOFs','energy_grids_ASCII')
mof_path = os.path.join(grid_path,'AHOKIR01-O.cif') #path to CIF of MOF

mof = read(mof_path)
site_idx = [atom.index for atom in mof if atom.symbol == 'O'][-1]
ads = adsorbate_constructor(ads='CH4',d_MX1=3.0)
mof_adsorbate = ads.get_adsorbate_grid(atoms_path=mof_path,
	site_idx=site_idx,grid_path=grid_path,grid_format='ASCII')