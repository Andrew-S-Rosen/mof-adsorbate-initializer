import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Fe-MOF-74.cif') #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_length = 2.0 #desired distance between site_index and ads_species
NN_method = 'crystal' #Pymatgen algorithm to detect bonding environment

#add adsorbate
ads_species = 'HOH' #adsorbate species (eta1-O)
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx)
new_mof_atoms, new_mof_name = ads.get_adsorbate_pm(mof_path,NN_method=NN_method,connect=2)
