import os
from mai.adsorbate_constructor import adsorbate_constructor
from ase.io import read

mof_path = os.path.join('example_MOFs','AHOKIR01-O.cif') #path to CIF of MOF
grid_path = os.path.join('example_MOFs','energy_grids_ASCII')
adsorbate = 'CH4' #adsorbate species
site = 'O' #adsorption site
max_dist = 3.0 #maximum distance between O and C

mof = read(mof_path)
site_idx = [atom.index for atom in mof if atom.symbol == site][-1]
ads = adsorbate_constructor(adsorbate,max_dist,site_idx=site_idx)
mof_adsorbate = ads.get_adsorbate_grid(mof_path,grid_path=grid_path,grid_format='ASCII')