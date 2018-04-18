import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'examples/oxygenated_MOFs/'
new_mofs_path = 'examples/add_H/'
site_species = 'O'
ads_species = 'H'
bond_length = 1.0
NN_method = 'vire'

for filename in os.listdir(mof_path):
	filepath = os.path.join(mof_path,filename)
	ads_const = adsorbate_constructor(ads_species,bond_length,
		site_species=site_species)
	mof_adsorbate, mof_name = ads_const.get_adsorbate_pm(filepath,NN_method,
		new_mofs_path=new_mofs_path)