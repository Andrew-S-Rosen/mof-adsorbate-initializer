from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Cu-BTC.cif' #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_dist = 1.75 #desired distance between site_idx and ads_species
ads_species = 'O' #adsorbate species

#add adsorbate
ads = adsorbate_constructor(ads_species,bond_dist,site_idx=site_idx)
new_mof_atoms = ads.get_adsorbate(mof_path)