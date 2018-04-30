import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'oxygenated_MOFs/'
new_mofs_path = 'add_N2/'
max_dist = 3.5
overlap_tol = 1.3
mol_species = 'N2'
site_species = 'O'
for filename in os.listdir(mof_path):
	filepath = os.path.join(mof_path,filename)
	ads_const = adsorbate_constructor(mol_species,max_dist,
		site_species=site_species,overlap_tol=overlap_tol)
	mof_adsorbate, mof_name = ads_const.get_adsorbate_raspa(filepath,
		new_mofs_path=new_mofs_path)
