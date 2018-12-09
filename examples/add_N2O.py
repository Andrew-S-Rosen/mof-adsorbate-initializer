import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','ANUGIA.cif') #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_length = 2.0 #desired distance between site_index and ads_species
neighbor_algo = 'okeeffe' #Pymatgen algorithm to detect bonding environment

#add adsorbate
ads_species = 'ON2' #adsorbate species (eta1-O)
NN_dist = 1.13
NO_dist = 1.19
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx)
new_mof_atoms, new_mof_name = ads.get_adsorbate_pm(mof_path,NN_method=neighbor_algo,
	d_bond=NO_dist,d_bond2=NN_dist)

ads_species = 'N2O' #adsorbate species (eta1-N)
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx)
new_mof_atoms, new_mof_name = ads.get_adsorbate_pm(mof_path,NN_method=neighbor_algo,
	d_bond=NN_dist,d_bond2=NO_dist)
