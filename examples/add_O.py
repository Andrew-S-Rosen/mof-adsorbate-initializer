import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'bare_MOFs/'
new_mofs_path = 'add_O/'
oms_data_path = os.path.join(mof_path,'OMS_data')
ads_species = 'O'
bond_length = 2.0
overlap_tol = 1.3

for filename in os.listdir(mof_path):
	filepath = os.path.join(mof_path,filename)
	ads_const = adsorbate_constructor(ads_species,bond_length,
		overlap_tol=overlap_tol)
	mof_adsorbate_list, mof_name_list = ads_const.get_adsorbate_zeo_oms(filepath,
		oms_data_path=oms_data_path,new_mofs_path=new_mofs_path)