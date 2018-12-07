import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','ANUGIA.cif') #path to CIF of MOF
site_idx = 0 #index of adsorption site
ads_species = 'O2' #adsorbate species
bond_length = 2.0 #desired distance between site_index and ads_species
NN_method = 'okeeffe' #Pymatgen algorithm to detect bonding environment

#add adsorbate
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx)
new_mof_atoms, new_mof_name = ads.get_adsorbate_pm(mof_path,NN_method=NN_method,eta=1)
