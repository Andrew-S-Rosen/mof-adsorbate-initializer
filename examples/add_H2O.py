import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Ni-BTP.cif') #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_length = 2.0 #desired distance between site_index and ads_species
OH_dist = 0.96

#add adsorbate
ads_species = 'HOH' #adsorbate species (eta1-O)
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx,
	d_bond=OH_dist,d_bond2=OH_dist,angle=120,angle2=104.5,connect=2)
new_mof_atoms, new_mof_name = ads.get_adsorbate_pm(mof_path)
