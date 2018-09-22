import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','ANUGIA.cif') #path to CIF of MOF
oms_path = os.path.join('example_MOFs','oms_results') #path to .oms(ex) data
ads_species = 'O' #adsorbate species
bond_length = 2.0 #desired distance between site_index and ads_species

#add adsorbate
ads = adsorbate_constructor(ads_species,bond_length)
new_mof_atoms, new_mof_name = ads.get_adsorbate_oms(mof_path,
	oms_data_path=oms_path,oms_format='omd')