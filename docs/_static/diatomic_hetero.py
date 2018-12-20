from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Ni-BTP.cif' #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_length = 1.5 #desired distance between site_index and ads_species
CO_length = 1.13 #C-O bond length

#add adsorbate
ads_species = 'CO' #adsorbate species (eta1-C)
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx,
	d_bond=CO_length)
new_mof_atoms1 = ads.get_adsorbate(mof_path)

ads_species = 'OC' #adsorbate species (eta1-O)
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx,
	d_bond=CO_length)
new_mof_atoms2 = ads.get_adsorbate(mof_path)
