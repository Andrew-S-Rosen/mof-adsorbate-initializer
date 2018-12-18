import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Ni-BTP.cif') #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_length = 1.75 #desired distance between site_index and ads_species

#add adsorbate
ads_species = 'O2_end' #adsorbate species
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx,eta=1)
new_mof_atoms = ads.get_adsorbate(mof_path)

ads_species = 'O2_side' #adsorbate species
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx,eta=2)
new_mof_atoms = ads.get_adsorbate(mof_path)
