from mai.adsorbate_constructor import adsorbate_constructor
from ase.io import read

grid_path = 'energy_grids' #path to energy grids
mof_path = 'mymof.cif' #path to CIF of MOF
site_idx = -1 #adsorption site

mof = read(mof_path)
ads = adsorbate_constructor(ads='CH4',d_MX1=3.0,site_idx=site_idx)
mof_adsorbate = ads.get_adsorbate_grid(atoms_path=mof_path,
	grid_path=grid_path,grid_format='ASCII')