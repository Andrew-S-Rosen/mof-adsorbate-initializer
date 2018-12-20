from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Ni2Cl2-BTDD.cif' #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_dist = 2.0 #desired distance between site_index and ads_species
NN_dist = 1.13 #N-N distance in N2O
NO_dist = 1.19 #N-O distance in N2O

#add adsorbate
ads_species = 'N2O' #adsorbate species (eta1-N)
ads = adsorbate_constructor(ads_species,bond_dist,site_idx=site_idx,
	d_bond=NN_dist,d_bond2=NO_dist)
new_mof_atoms1 = ads.get_adsorbate(mof_path)

ads_species = 'ON2' #adsorbate species (eta1-O)
ads = adsorbate_constructor(ads_species,bond_dist,site_idx=site_idx,
	d_bond=NO_dist,d_bond2=NN_dist)
new_mof_atoms2 = ads.get_adsorbate(mof_path)
