from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Ni-BTP.cif' #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_dist = 1.5 #desired distance between site_index and ads_species
OO_dist = 1.2 #O-O bond distance in O2
end_on_angle = 120 #angle for end-on adsorption
side_on_angle = 90 #angle for side-on adsorption

#add adsorbate
ads_species = 'O2_end' #adsorbate species
ads = adsorbate_constructor(ads_species,bond_dist,site_idx=site_idx,
	eta=1,d_bond=OO_dist,angle=end_on_angle)
new_mof_atoms1 = ads.get_adsorbate(mof_path)

ads_species = 'O2_side' #adsorbate species
ads = adsorbate_constructor(ads_species,bond_dist,site_idx=site_idx,
	eta=2,d_bond=OO_dist,angle=side_on_angle)
new_mof_atoms2 = ads.get_adsorbate(mof_path)
