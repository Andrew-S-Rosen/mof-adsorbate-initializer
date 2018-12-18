import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Mn-DSBDC.cif') #path to CIF of MOF
site_idx = 6 #index of adsorption site
bond_length = 1.75 #desired distance between site_index and ads_species
ads_species = 'O' #adsorbate species

#add adsorbate
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx)
new_mof_atoms = ads.get_adsorbate(mof_path)
