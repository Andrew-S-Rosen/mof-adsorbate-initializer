from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Sc-MIL-88B.cif' #path to CIF of MOF
site_idx = 0 #index of adsorption site
bond_length = 2.0 #desired distance between site_index and ads_species
ads_species = 'HOH' #adsorbate species
OH_dist = 0.96 #O-H distance in HOH
angle1 = 120 #M-O-H angle
angle2 = 104.5 #H-O-H angle
connect_atom = 2 #adsorption via O atom in HOH

#add adsorbate
ads = adsorbate_constructor(ads_species,bond_length,site_idx=site_idx,
	d_bond=OH_dist,d_bond2=OH_dist,angle=angle1,angle2=angle2,
	connect=connect_atom)
new_mof_atoms = ads.get_adsorbate(mof_path)
